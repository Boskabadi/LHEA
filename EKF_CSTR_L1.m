function EKF_CSTR_L1(block)
% L1 EKF (no ML): uses measurements [V, CX, DO] fast + MEV assay every 12 h
% InputPort(1): x_true (6)  (unused, kept for compatibility)
% InputPort(2): y_meas (4)  = [V; CX; DO; MEV_assay]
% InputPort(3): dt (1)      = hours
    setup(block);
end

function setup(block)
    block.NumInputPorts  = 3;
    block.NumOutputPorts = 1;

    block.InputPort(1).Dimensions = 6;   % x_true (unused)
    block.InputPort(2).Dimensions = 4;   % [V, CX, DO, MEV]
    block.InputPort(3).Dimensions = 1;   % dt [h]

    for i = 1:3
        block.InputPort(i).DatatypeID = 0;          % double
        block.InputPort(i).DirectFeedthrough = true;
    end

    block.OutputPort(1).Dimensions = 6;  % xhat
    block.OutputPort(1).DatatypeID = 0;

    block.SampleTimes = [0 0];
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('PostPropagationSetup', @PostPropSetup);
    block.RegBlockMethod('InitializeConditions', @InitConditions);
    block.RegBlockMethod('Outputs',              @Output);
    block.RegBlockMethod('Update',               @Update);
end

function PostPropSetup(block)
    % Dwork list:
    % 1 xhat(6)
    % 2 P(36)
    % 3 k
    % 4 n_assay
    % 5 sum_abs_err
    % 6 sum_sq_err
    % 7 mean_y (for R2)
    % 8 M2_y   (for R2)  SST = M2_y
    % 9 nis_in
    % 10 nis_total
    block.NumDworks = 10;

    block.Dwork(1).Name='xhat';
    block.Dwork(1).Dimensions=6;
    block.Dwork(1).DatatypeID=0;
    block.Dwork(1).Complexity='Real';
    block.Dwork(1).UsedAsDiscState=true;

    block.Dwork(2).Name='P';
    block.Dwork(2).Dimensions=36;
    block.Dwork(2).DatatypeID=0;
    block.Dwork(2).Complexity='Real';
    block.Dwork(2).UsedAsDiscState=true;

    block.Dwork(3).Name='k';
    block.Dwork(3).Dimensions=1;
    block.Dwork(3).DatatypeID=0;
    block.Dwork(3).Complexity='Real';
    block.Dwork(3).UsedAsDiscState=true;

    block.Dwork(4).Name='n_assay';
    block.Dwork(4).Dimensions=1;
    block.Dwork(4).DatatypeID=0;
    block.Dwork(4).Complexity='Real';
    block.Dwork(4).UsedAsDiscState=true;

    block.Dwork(5).Name='sum_abs_err';
    block.Dwork(5).Dimensions=1;
    block.Dwork(5).DatatypeID=0;
    block.Dwork(5).Complexity='Real';
    block.Dwork(5).UsedAsDiscState=true;

    block.Dwork(6).Name='sum_sq_err';
    block.Dwork(6).Dimensions=1;
    block.Dwork(6).DatatypeID=0;
    block.Dwork(6).Complexity='Real';
    block.Dwork(6).UsedAsDiscState=true;

    block.Dwork(7).Name='mean_y';
    block.Dwork(7).Dimensions=1;
    block.Dwork(7).DatatypeID=0;
    block.Dwork(7).Complexity='Real';
    block.Dwork(7).UsedAsDiscState=true;

    block.Dwork(8).Name='M2_y';
    block.Dwork(8).Dimensions=1;
    block.Dwork(8).DatatypeID=0;
    block.Dwork(8).Complexity='Real';
    block.Dwork(8).UsedAsDiscState=true;

    block.Dwork(9).Name='nis_in';
    block.Dwork(9).Dimensions=1;
    block.Dwork(9).DatatypeID=0;
    block.Dwork(9).Complexity='Real';
    block.Dwork(9).UsedAsDiscState=true;

    block.Dwork(10).Name='nis_total';
    block.Dwork(10).Dimensions=1;
    block.Dwork(10).DatatypeID=0;
    block.Dwork(10).Complexity='Real';
    block.Dwork(10).UsedAsDiscState=true;
end

function InitConditions(block)
    x0 = [5000; 106; 14.8; 2.52; 1.2; 0.008];
    block.Dwork(1).Data = x0;

    P0 = diag([1e2 5 5 5 5 1e-4]);
    block.Dwork(2).Data = reshape(P0,36,1);

    block.Dwork(3).Data  = 0;

    % metrics init
    block.Dwork(4).Data  = 0;    % n_assay
    block.Dwork(5).Data  = 0;    % sum_abs_err
    block.Dwork(6).Data  = 0;    % sum_sq_err
    block.Dwork(7).Data  = 0;    % mean_y
    block.Dwork(8).Data  = 0;    % M2_y
    block.Dwork(9).Data  = 0;    % nis_in
    block.Dwork(10).Data = 0;    % nis_total
end

function Output(block)
    block.OutputPort(1).Data = block.Dwork(1).Data;
end

function Update(block)
    %#ok<NASGU>
    x_true = block.InputPort(1).Data;      % not used
    y      = block.InputPort(2).Data;      % [V; CX; DO; MEV_assay]
    dt     = block.InputPort(3).Data;      % hours

    xhat = block.Dwork(1).Data;
    P    = reshape(block.Dwork(2).Data,6,6);
    k    = block.Dwork(3).Data;

    % ---- metrics states ----
    n_assay    = block.Dwork(4).Data;
    sum_abs_e  = block.Dwork(5).Data;
    sum_sq_e   = block.Dwork(6).Data;
    mean_y     = block.Dwork(7).Data;
    M2_y       = block.Dwork(8).Data;
    nis_in     = block.Dwork(9).Data;
    nis_total  = block.Dwork(10).Data;

    % ---------- parameters for CSTRO_ode + linearize_CSTRO ----------
    d = 0.9;
    p.A_xsec = pi/4*d^2;
    p.gamma0 = 0;
    p.gamma1 = 1;
    p.Y_O2_LAC_growth = 0.52;
    p.Y_O2_LAC_prod   = 0.12;

    % ---------- known feed ----------
    u = [88; 16; 3; 0.4; 60; 60];

    % ---------- noises ----------
    Q = diag([1e-2 1e-3 1e-3 1e-4 1e-4 1e-6]);

    % Measurement noise
    sigma_V   = 10;
    sigma_CX  = 0.5;
    sigma_DO  = 5e-4;
    sigma_MEV = 0.05;   % (g/L) tune

    R_V   = sigma_V^2;
    R_CX  = sigma_CX^2;
    R_DO  = sigma_DO^2;
    R_MEV = sigma_MEV^2;

    % ---------- measurement models ----------
    H_V   = [1 0 0 0 0 0];
    H_CX  = [0 1 0 0 0 0];
    H_DO  = [0 0 0 0 0 1];
    H_MEV = [0 0 0 0 1 0];

    % ---------- cadences ----------
    dt_sec = dt*3600;
    Ns_Level = max(1, round(5          / dt_sec));
    Ns_CX    = max(1, round(10         / dt_sec));
    Ns_DO    = max(1, round(5          / dt_sec));
    Ns_MEV   = max(1, round((12*3600)  / dt_sec)); % 12 h

    % ---------- prediction ----------
    f = CSTRO_ode(xhat, u, p);
    xhat = xhat + dt*f;

    [A,~,~,~,~,~] = linearize_CSTRO(xhat, u, p);
    Ad = expm(A*dt);
    P  = Ad*P*Ad' + Q;

    % Save predicted copies (for MEV NIS)
    xhat_pred = xhat;
    P_pred    = P;

    % ---------- multi-rate update ----------
    Hk = [];
    yk = [];
    Rk = [];
    mev_included = false;

    if mod(k, Ns_Level)==0
        Hk = [Hk; H_V];  yk = [yk; y(1)];  Rk = blkdiag(Rk, R_V);
    end
    if mod(k, Ns_CX)==0
        Hk = [Hk; H_CX]; yk = [yk; y(2)];  Rk = blkdiag(Rk, R_CX);
    end
    if mod(k, Ns_DO)==0
        Hk = [Hk; H_DO]; yk = [yk; y(3)];  Rk = blkdiag(Rk, R_DO);
    end
    if mod(k, Ns_MEV)==0
        Hk = [Hk; H_MEV]; yk = [yk; y(4)]; Rk = blkdiag(Rk, R_MEV);
        mev_included = true;
    end

    if ~isempty(Hk)
        innov = yk - Hk*xhat;
        S = Hk*P*Hk' + Rk;
        K = P*Hk'/S;

        xhat = xhat + K*innov;
        P    = (eye(6) - K*Hk)*P;
    end

    % ---------- metrics update ONLY when MEV assay exists ----------
    if mev_included
        % NIS for MEV channel (scalar, using predicted covariance)
        innov_mev = y(4) - (H_MEV*xhat_pred);
        S_mev     = H_MEV*P_pred*H_MEV' + R_MEV;
        nis       = (innov_mev*innov_mev) / S_mev;

        % 95% bounds for chi-square(dof=1): [0.025, 0.975]
        nis_lo = 0.000982;   % chi2inv(0.025,1)
        nis_hi = 5.023886;   % chi2inv(0.975,1)

        nis_total = nis_total + 1;
        if (nis >= nis_lo) && (nis <= nis_hi)
            nis_in = nis_in + 1;
        end

        % error AFTER update (what you care about)
        e = y(4) - xhat(5);

        n_assay   = n_assay + 1;
        sum_abs_e = sum_abs_e + abs(e);
        sum_sq_e  = sum_sq_e  + e*e;

        % Welford update for assay variance (SST)
        if n_assay == 1
            mean_y = y(4);
            M2_y   = 0;
        else
            dy    = y(4) - mean_y;
            mean_y = mean_y + dy / n_assay;
            dy2   = y(4) - mean_y;
            M2_y  = M2_y + dy*dy2;
        end

        % compute metrics safely
        rmse = sqrt(sum_sq_e / max(1,n_assay));
        mae  = sum_abs_e / max(1,n_assay);
        sst  = M2_y;                   % sum (y - mean)^2 over assays
        sse  = sum_sq_e;
        if sst > 0
            R2 = 1 - (sse / sst);
        else
            R2 = NaN;
        end
        nis_cov = nis_in / max(1,nis_total);

        assignin('base','L1_metrics',[rmse, mae, R2, nis_cov]);
        assignin('base','L1_nAssay',n_assay);
    end

    % ---------- save back ----------
    block.Dwork(1).Data = xhat;
    block.Dwork(2).Data = reshape(P,36,1);
    block.Dwork(3).Data = k + 1;

    block.Dwork(4).Data  = n_assay;
    block.Dwork(5).Data  = sum_abs_e;
    block.Dwork(6).Data  = sum_sq_e;
    block.Dwork(7).Data  = mean_y;
    block.Dwork(8).Data  = M2_y;
    block.Dwork(9).Data  = nis_in;
    block.Dwork(10).Data = nis_total;

    % --- Metrics on LOV using x_true (plant truth) ---
    lov_true = x_true(5);
    lov_hat  = xhat(5);
    e = lov_true - lov_hat;

    persistent SSE SAE N
    if isempty(SSE), SSE=0; SAE=0; N=0; end

    SSE = SSE + e^2;
    SAE = SAE + abs(e);
    N   = N + 1;

    RMSE = sqrt(SSE / max(N,1));
    MAE  = SAE / max(N,1);

    % R2: only if truth varies
    % (if truth is constant, R2 is undefined -> keep NaN)
    persistent sumY sumY2
    if isempty(sumY), sumY=0; sumY2=0; end
    sumY  = sumY  + lov_true;
    sumY2 = sumY2 + lov_true^2;

    meanY = sumY / max(N,1);
    SST = sumY2 - N*meanY^2;

    if SST > 1e-12
        R2 = 1 - (SSE / SST);
    else
        R2 = NaN;  % constant truth => undefined
    end

    assignin('base','L1_metrics',[RMSE MAE R2]);

    end
