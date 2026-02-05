function EKF_CSTR_L3(block)
    setup(block);
end

function setup(block)
    block.NumInputPorts  = 3;
    block.NumOutputPorts = 1;

    block.InputPort(1).Dimensions = 6;   % x_true (SIM only)
    block.InputPort(2).Dimensions = 4;   % [V, CX, DO, MEV_assay]
    block.InputPort(3).Dimensions = 1;   % dt [h]

    for i = 1:3
        block.InputPort(i).DatatypeID = 0; % double
        block.InputPort(i).DirectFeedthrough = true;
    end

    block.OutputPort(1).Dimensions = 7;  % 6 states + qscale
    block.OutputPort(1).DatatypeID = 0;

    block.SampleTimes = [-1 0];  % inherit
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('PostPropagationSetup', @PostPropSetup);
    block.RegBlockMethod('InitializeConditions', @InitConditions);
    block.RegBlockMethod('Outputs',              @Output);
    block.RegBlockMethod('Update',               @Update);
end

function PostPropSetup(block)
    block.NumDworks = 8;

    % xhat
    block.Dwork(1).Name='xhat';
    block.Dwork(1).Dimensions=6;
    block.Dwork(1).DatatypeID=0;
    block.Dwork(1).Complexity='Real';
    block.Dwork(1).UsedAsDiscState=true;

    % P
    block.Dwork(2).Name='P';
    block.Dwork(2).Dimensions=36;
    block.Dwork(2).DatatypeID=0;
    block.Dwork(2).Complexity='Real';
    block.Dwork(2).UsedAsDiscState=true;

    % k
    block.Dwork(3).Name='k';
    block.Dwork(3).Dimensions=1;
    block.Dwork(3).DatatypeID=0;
    block.Dwork(3).Complexity='Real';
    block.Dwork(3).UsedAsDiscState=true;

    % qscale
    block.Dwork(4).Name='qscale';
    block.Dwork(4).Dimensions=1;
    block.Dwork(4).DatatypeID=0;
    block.Dwork(4).Complexity='Real';
    block.Dwork(4).UsedAsDiscState=true;

    % tA (kept for compatibility)
    block.Dwork(5).Name='tA';
    block.Dwork(5).Dimensions=1;
    block.Dwork(5).DatatypeID=0;
    block.Dwork(5).Complexity='Real';
    block.Dwork(5).UsedAsDiscState=true;

    % last MEV assay
    block.Dwork(6).Name='mev_last';
    block.Dwork(6).Dimensions=1;
    block.Dwork(6).DatatypeID=0;
    block.Dwork(6).Complexity='Real';
    block.Dwork(6).UsedAsDiscState=true;

    % theta (3x1)
    block.Dwork(7).Name='theta';
    block.Dwork(7).Dimensions=3;
    block.Dwork(7).DatatypeID=0;
    block.Dwork(7).Complexity='Real';
    block.Dwork(7).UsedAsDiscState=true;

    % Ptheta (3x3 flattened)
    block.Dwork(8).Name='Ptheta';
    block.Dwork(8).Dimensions=9;
    block.Dwork(8).DatatypeID=0;
    block.Dwork(8).Complexity='Real';
    block.Dwork(8).UsedAsDiscState=true;
end

function InitConditions(block)
    block.Dwork(1).Data = [5000; 106; 14.8; 2.52; 1.2; 0.008];
    block.Dwork(2).Data = reshape(diag([1e2 5 5 5 5 1e-4]),36,1);
    block.Dwork(3).Data = 0;
    block.Dwork(4).Data = 1.0;
    block.Dwork(5).Data = 0.0;
    block.Dwork(6).Data = 1.2;

    block.Dwork(7).Data = [0;0;0];
    block.Dwork(8).Data = reshape(1e4*eye(3),9,1);

    % reset logs in base workspace
    assignin('base','L3_time_hr',[]);
    assignin('base','L3_mev_true',[]);
    assignin('base','L3_mev_hat',[]);
    assignin('base','L3_metrics',[NaN NaN NaN NaN]);
end

function Output(block)
    xhat = block.Dwork(1).Data;
    qscale = block.Dwork(4).Data;
    block.OutputPort(1).Data = [xhat; qscale];
end

function Update(block)
    x_true = block.InputPort(1).Data;   % SIM only
    y      = block.InputPort(2).Data;   % [V, CX, DO, MEV_assay]
    dt     = block.InputPort(3).Data;   % [h]

    xhat   = block.Dwork(1).Data;
    P      = reshape(block.Dwork(2).Data,6,6);
    k      = block.Dwork(3).Data;
    qscale = block.Dwork(4).Data;
    tA     = block.Dwork(5).Data; %#ok<NASGU>
    mev_last = block.Dwork(6).Data;

    theta  = block.Dwork(7).Data;           % 3x1
    Ptheta = reshape(block.Dwork(8).Data,3,3);

    if ~(isfinite(dt) && dt > 0)
        dt = 1/3600;
    end

    % ---- params required by linearize_CSTRO ----
    d = 0.9;
    p.A_xsec = pi/4*d^2;
    p.gamma0 = 0; p.gamma1 = 1;
    p.Y_O2_LAC_growth = 0.52;
    p.Y_O2_LAC_prod   = 0.12;

    % ---- known feed ----
    u = [88; 16; 3.0; 0.4; 60; 60];

    % ---- noises ----
    Q = diag([1e-2 1e-3 1e-3 1e-4 1e-4 1e-6]);

    sigV=10; sigCX=0.5; sigDO=5e-4; sigMEV=0.05;

    % ---- fast sensor rates ----
    dt_sec = dt*3600;
    Ns_V  = max(1, round(5  / dt_sec));
    Ns_CX = max(1, round(10 / dt_sec));
    Ns_DO = max(1, round(5  / dt_sec));

    % ---- sparse assay schedule (12h) ----
    T_assay_hr = 12;
    t_hr = block.CurrentTime;
    assay_now = (mod(t_hr, T_assay_hr) < dt);

    % ---------- prediction ----------
    p.q_max_scale = qscale;
    f = CSTRO_ode(xhat, u, p);

    % ML residual on MEV rate: r = theta' * [1; CX; DO]
    phi = [1; xhat(2); xhat(6)];
    r_mev = theta.' * phi;                  % [g/L/h]
    r_mev = max(min(r_mev, 5e-4), -5e-4);   % bound
    f(5) = f(5) + r_mev;

    xhat_pred = xhat + dt*f;   % store predicted for learning/stat
    xhat = xhat_pred;

    [A,~,~,~,~,~] = linearize_CSTRO(xhat, u, p);
    Ad = expm(A*dt);
    P  = Ad*P*Ad' + Q;

    % ---------- measurement update ----------
    Hk=[]; Rk=[]; yk=[];

    if mod(k,Ns_V)==0
        Hk=[Hk; 1 0 0 0 0 0]; Rk=blkdiag(Rk,sigV^2);  yk=[yk; y(1)];
    end
    if mod(k,Ns_CX)==0
        Hk=[Hk; 0 1 0 0 0 0]; Rk=blkdiag(Rk,sigCX^2); yk=[yk; y(2)];
    end
    if mod(k,Ns_DO)==0
        Hk=[Hk; 0 0 0 0 0 1]; Rk=blkdiag(Rk,sigDO^2); yk=[yk; y(3)];
    end
    if assay_now
        Hk=[Hk; 0 0 0 0 1 0]; Rk=blkdiag(Rk,sigMEV^2); yk=[yk; y(4)];
        mev_last = y(4);
    end

    if ~isempty(Hk)
        innov = yk - Hk*xhat;
        S = Hk*P*Hk' + Rk;
        K = P*Hk'/S;
        xhat = xhat + K*innov;
        P = (eye(6)-K*Hk)*P;
    end

    % ---------- L3 learning (RLS) only at assay ----------
    if assay_now
        % learning signal = assay - predicted (pre-update)
        e_mev = mev_last - xhat_pred(5);

        % RLS
        lambda = 0.995;
        denom = lambda + (phi.'*Ptheta*phi);
        Kth = (Ptheta*phi)/denom;
        theta = theta + Kth * e_mev;
        Ptheta = (Ptheta - Kth*(phi.'*Ptheta))/lambda;

        % qscale adaptation (slow, bounded)
        eta_q = 2e-3;
        qscale = qscale + eta_q * e_mev;
        qscale = min(max(qscale,0.7),1.3);

        % ---------- STATISTICS (MEV only, at assay times) ----------
        tlog = evalin('base','L3_time_hr');
        yT   = evalin('base','L3_mev_true');
        yH   = evalin('base','L3_mev_hat');

        tlog = [tlog; t_hr];
        yT   = [yT;   x_true(5)];
        yH   = [yH;   xhat(5)];

        e = yT - yH;
        rmse = sqrt(mean(e.^2));
        mae  = mean(abs(e));
        bias = mean(e);
        maxe = max(abs(e));

        assignin('base','L3_time_hr',tlog);
        assignin('base','L3_mev_true',yT);
        assignin('base','L3_mev_hat',yH);
        assignin('base','L3_metrics',[rmse mae bias maxe]);
    end

    % save back
    block.Dwork(1).Data = xhat;
    block.Dwork(2).Data = reshape(P,36,1);
    block.Dwork(3).Data = k+1;
    block.Dwork(4).Data = qscale;
    block.Dwork(5).Data = 0.0;      % kept
    block.Dwork(6).Data = mev_last;
    block.Dwork(7).Data = theta;
    block.Dwork(8).Data = reshape(Ptheta,9,1);
end
