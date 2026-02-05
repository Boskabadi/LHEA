function EKF_CSTR_L2_S9(block)
% L2 EKF (bias-robust): Scenario 9
% Measurements: [V, CX, DO] fast + MEV assay every 12 h
% Augmented state: [V CX LAC N MEV DO  bV  bDO]
%
% InputPort(1): x_true (6)   (used ONLY for evaluation in simulation)
% InputPort(2): y_meas (4)   = [V; CX; DO; MEV_assay]
% InputPort(3): dt (1)       = hours
%
% OutputPort(1): xhat_aug (8)

    setup(block);
end

function setup(block)
    block.NumInputPorts  = 3;
    block.NumOutputPorts = 1;

    block.InputPort(1).Dimensions = 6;   % x_true (for metrics only)
    block.InputPort(2).Dimensions = 4;   % [V CX DO MEV]
    block.InputPort(3).Dimensions = 1;   % dt [h]

    for i = 1:3
        block.InputPort(i).DatatypeID = 0;          % double
        block.InputPort(i).DirectFeedthrough = true;
    end

    block.OutputPort(1).Dimensions = 8;             % augmented state
    block.OutputPort(1).DatatypeID = 0;

    block.SampleTimes = [0 0];
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('PostPropagationSetup', @PostPropSetup);
    block.RegBlockMethod('InitializeConditions', @InitConditions);
    block.RegBlockMethod('Outputs',              @Output);
    block.RegBlockMethod('Update',               @Update);
end

function PostPropSetup(block)
    block.NumDworks = 3;

    % xhat_aug
    block.Dwork(1).Name = 'xhat';
    block.Dwork(1).Dimensions = 8;
    block.Dwork(1).DatatypeID = 0;
    block.Dwork(1).Complexity = 'Real';
    block.Dwork(1).UsedAsDiscState = true;

    % P
    block.Dwork(2).Name = 'P';
    block.Dwork(2).Dimensions = 64;
    block.Dwork(2).DatatypeID = 0;
    block.Dwork(2).Complexity = 'Real';
    block.Dwork(2).UsedAsDiscState = true;

    % counter k
    block.Dwork(3).Name = 'k';
    block.Dwork(3).Dimensions = 1;
    block.Dwork(3).DatatypeID = 0;
    block.Dwork(3).Complexity = 'Real';
    block.Dwork(3).UsedAsDiscState = true;
end

function InitConditions(block)
    x_phys0 = [5000; 106; 14.8; 2.52; 1.2; 0.008];
    b0      = [0; 0];  % [bV; bDO]
    block.Dwork(1).Data = [x_phys0; b0];

    P_phys = diag([1e2 5 5 5 5 1e-4]);
    P_bias = diag([1e4 1e4]);     % allow bias to move
    P0 = blkdiag(P_phys, P_bias);

    block.Dwork(2).Data = reshape(P0,64,1);
    block.Dwork(3).Data = 0;

    % reset metrics in base workspace
    assignin('base','L2_metrics',[NaN NaN NaN NaN]); % [RMSE MAE R2 MaxAbs]
end

function Output(block)
    block.OutputPort(1).Data = block.Dwork(1).Data;
end

function Update(block)
    x_true = block.InputPort(1).Data;   % evaluation only (SIM)
    y      = block.InputPort(2).Data;   % [V; CX; DO; MEV_assay]
    dt     = block.InputPort(3).Data;   % [h]

    xhat_a = block.Dwork(1).Data;       % 8x1
    P      = reshape(block.Dwork(2).Data,8,8);
    k      = block.Dwork(3).Data;

    % split
    x_phys = xhat_a(1:6);
    bV     = xhat_a(7);
    bDO    = xhat_a(8);

    % ---- parameters (must match linearize_CSTRO expectations) ----
    d = 0.9;
    p.A_xsec = pi/4*d^2;
    p.gamma0 = 0;
    p.gamma1 = 1;
    p.Y_O2_LAC_growth = 0.52;
    p.Y_O2_LAC_prod   = 0.12;

    % ---- known feed ----
    u = [88; 16; 3; 0.4; 60; 60];

    % ---- process noise ----
    Q_phys = diag([1e-2 1e-3 1e-3 1e-4 1e-4 1e-6]);
    Q_bias = diag([1e-3 1e-3]);
    Q_a    = blkdiag(Q_phys, Q_bias);

    % ---- measurement noise ----
    sigma_V   = 10;
    sigma_CX  = 0.5;
    sigma_DO  = 5e-4;
    sigma_MEV = 0.05;

    R_V   = sigma_V^2;
    R_CX  = sigma_CX^2;
    R_DO  = sigma_DO^2;
    R_MEV = sigma_MEV^2;

    % ---- measurement rows (with biases) ----
    H_V   = [1 0 0 0 0 0 1 0];
    H_CX  = [0 1 0 0 0 0 0 0];
    H_DO  = [0 0 0 0 0 1 0 1];
    H_MEV = [0 0 0 0 1 0 0 0];

    % ---- cadences ----
    dt_sec   = dt*3600;
    Ns_Level = max(1, round(5 / dt_sec));
    Ns_CX    = max(1, round(10 / dt_sec));
    Ns_DO    = max(1, round(5 / dt_sec));
    Ns_MEV   = max(1, round((12*3600) / dt_sec));

    % ---------- prediction ----------
    f = CSTRO_ode(x_phys, u, p);
    x_phys = x_phys + dt*f;

    [A,~,~,~,~,~] = linearize_CSTRO(x_phys, u, p);

    A_a = zeros(8,8);
    A_a(1:6,1:6) = A;   % biases random-walk

    Ad_a = expm(A_a*dt);
    P    = Ad_a*P*Ad_a' + Q_a;

    xhat_a = [x_phys; bV; bDO];

    % ---------- multi-rate measurement update ----------
    Hk = [];
    yk = [];
    Rk = [];

    if mod(k, Ns_Level)==0
        Hk = [Hk; H_V];
        yk = [yk; y(1)];
        Rk = blkdiag(Rk, R_V);
    end
    if mod(k, Ns_CX)==0
        Hk = [Hk; H_CX];
        yk = [yk; y(2)];
        Rk = blkdiag(Rk, R_CX);
    end
    if mod(k, Ns_DO)==0
        Hk = [Hk; H_DO];
        yk = [yk; y(3)];
        Rk = blkdiag(Rk, R_DO);
    end
    if mod(k, Ns_MEV)==0
        Hk = [Hk; H_MEV];
        yk = [yk; y(4)];
        Rk = blkdiag(Rk, R_MEV);
    end

    if ~isempty(Hk)
        innov = yk - Hk*xhat_a;
        S     = Hk*P*Hk' + Rk;
        K     = P*Hk'/S;

        xhat_a = xhat_a + K*innov;
        P      = (eye(8) - K*Hk)*P;
    end

    % ---------- metrics on MEV (state 5) vs x_true (SIM only) ----------
    % same style as L1: accumulate error statistics online
    persistent n sumAbs sumSq sumY sumY2 sumYYhat maxAbs
    if isempty(n)
        n = 0; sumAbs = 0; sumSq = 0;
        sumY = 0; sumY2 = 0; sumYYhat = 0;
        maxAbs = 0;
    end

    y_true = x_true(5);
    y_hat  = xhat_a(5);
    e      = y_true - y_hat;

    n = n + 1;
    sumAbs = sumAbs + abs(e);
    sumSq  = sumSq  + e^2;
    sumY   = sumY   + y_true;
    sumY2  = sumY2  + y_true^2;
    sumYYhat = sumYYhat + (y_true - y_hat); % (not used directly, kept simple)
    maxAbs = max(maxAbs, abs(e));

    RMSE = sqrt(sumSq / max(n,1));
    MAE  = sumAbs / max(n,1);

    % R2 vs TRUE trajectory (not measurement)
    % R2 = 1 - SSE/SST ; handle constant y_true (SST≈0)
    SSE = sumSq;
    meanY = sumY / max(n,1);
    SST = sumY2 - n*(meanY^2);
    if SST > 1e-12
        R2 = 1 - SSE/SST;
    else
        R2 = NaN;
    end

    assignin('base','L2_metrics',[RMSE MAE R2 maxAbs]);

    % save
    block.Dwork(1).Data = xhat_a;
    block.Dwork(2).Data = reshape(P,64,1);
    block.Dwork(3).Data = k + 1;
end
