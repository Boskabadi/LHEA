function EKF_CSTR_L11(block)
% L1 EKF (no ML): uses measurements [V, CX, DO] fast + MEV assay every 12 h
% InputPort(1): x_true (6)  (used only for trajectory metrics)
% InputPort(2): y_meas (4)  = [V; CX; DO; MEV_assay]
% InputPort(3): dt (1)      = hours
%
% IMPORTANT:
% This code assumes Simulink model time is also in HOURS.
% If your model time is in seconds, convert the sampling periods below.

    setup(block);
end

function setup(block)
    block.NumInputPorts  = 3;
    block.NumOutputPorts = 1;

    block.InputPort(1).Dimensions = 6;   % x_true
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
    % Dworks
    %
    %  1  xhat(6)
    %  2  P(36)
    %  3  next_V_time
    %  4  next_CX_time
    %  5  next_DO_time
    %  6  next_MEV_time
    %
    % assay metrics
    %  7  n_assay
    %  8  sum_abs_err_assay
    %  9  sum_sq_err_assay
    % 10  mean_y_assay
    % 11  M2_y_assay
    % 12  nis_in
    % 13  nis_total
    %
    % full trajectory metrics against x_true(5)
    % 14  traj_SSE
    % 15  traj_SAE
    % 16  traj_N
    % 17  traj_sumY
    % 18  traj_sumY2

    block.NumDworks = 18;

    for i = 1:18
        block.Dwork(i).DatatypeID = 0;
        block.Dwork(i).Complexity = 'Real';
        block.Dwork(i).UsedAsDiscState = true;
    end

    block.Dwork(1).Name = 'xhat';                block.Dwork(1).Dimensions = 6;
    block.Dwork(2).Name = 'P';                   block.Dwork(2).Dimensions = 36;
    block.Dwork(3).Name = 'next_V_time';         block.Dwork(3).Dimensions = 1;
    block.Dwork(4).Name = 'next_CX_time';        block.Dwork(4).Dimensions = 1;
    block.Dwork(5).Name = 'next_DO_time';        block.Dwork(5).Dimensions = 1;
    block.Dwork(6).Name = 'next_MEV_time';       block.Dwork(6).Dimensions = 1;

    block.Dwork(7).Name = 'n_assay';             block.Dwork(7).Dimensions = 1;
    block.Dwork(8).Name = 'sum_abs_err_assay';   block.Dwork(8).Dimensions = 1;
    block.Dwork(9).Name = 'sum_sq_err_assay';    block.Dwork(9).Dimensions = 1;
    block.Dwork(10).Name = 'mean_y_assay';       block.Dwork(10).Dimensions = 1;
    block.Dwork(11).Name = 'M2_y_assay';         block.Dwork(11).Dimensions = 1;
    block.Dwork(12).Name = 'nis_in';             block.Dwork(12).Dimensions = 1;
    block.Dwork(13).Name = 'nis_total';          block.Dwork(13).Dimensions = 1;

    block.Dwork(14).Name = 'traj_SSE';           block.Dwork(14).Dimensions = 1;
    block.Dwork(15).Name = 'traj_SAE';           block.Dwork(15).Dimensions = 1;
    block.Dwork(16).Name = 'traj_N';             block.Dwork(16).Dimensions = 1;
    block.Dwork(17).Name = 'traj_sumY';          block.Dwork(17).Dimensions = 1;
    block.Dwork(18).Name = 'traj_sumY2';         block.Dwork(18).Dimensions = 1;
end

function InitConditions(block)
    x0 = [5000; 106; 14.8; 2.52; 1.2; 0.008];
    block.Dwork(1).Data = x0;

    P0 = diag([1e2 5 5 5 5 1e-4]);
    block.Dwork(2).Data = reshape(P0, 36, 1);

    % Start measurement clocks at t = 0 so the first sample is accepted immediately
    block.Dwork(3).Data = 0;   % next_V_time
    block.Dwork(4).Data = 0;   % next_CX_time
    block.Dwork(5).Data = 0;   % next_DO_time
    block.Dwork(6).Data = 0;   % next_MEV_time

    % Assay metrics
    block.Dwork(7).Data  = 0;
    block.Dwork(8).Data  = 0;
    block.Dwork(9).Data  = 0;
    block.Dwork(10).Data = 0;
    block.Dwork(11).Data = 0;
    block.Dwork(12).Data = 0;
    block.Dwork(13).Data = 0;

    % Trajectory metrics
    block.Dwork(14).Data = 0;
    block.Dwork(15).Data = 0;
    block.Dwork(16).Data = 0;
    block.Dwork(17).Data = 0;
    block.Dwork(18).Data = 0;

    % Workspace diagnostics
    assignin('base', 'L1_assay_log', zeros(0,8));
    assignin('base', 'L1_assay_log_cols', ...
        {'t_h','y_mev','xpred_mev','xpost_mev','innov_mev','P55_pred','K55_mev','nis'});
    assignin('base', 'L1_assay_metrics', [NaN NaN NaN NaN]);  % [RMSE MAE R2 NIS_cov]
    assignin('base', 'L1_traj_metrics',  [NaN NaN NaN]);      % [RMSE MAE R2]
end

function Output(block)
    % Note:
    % In a Level-2 MATLAB S-function, Outputs is called before Update.
    % So the visible jump from the current measurement update may appear one step later.
    block.OutputPort(1).Data = block.Dwork(1).Data;
end

function Update(block)
    x_true = block.InputPort(1).Data;
    y      = block.InputPort(2).Data;      % [V; CX; DO; MEV_assay]
    dt     = block.InputPort(3).Data;      % hours
    dt     = max(dt, eps);

    t = block.CurrentTime;                 % assumed HOURS

    xhat = block.Dwork(1).Data;
    P    = reshape(block.Dwork(2).Data, 6, 6);

    next_V_time   = block.Dwork(3).Data;
    next_CX_time  = block.Dwork(4).Data;
    next_DO_time  = block.Dwork(5).Data;
    next_MEV_time = block.Dwork(6).Data;

    % Assay metrics
    n_assay       = block.Dwork(7).Data;
    sum_abs_e_a   = block.Dwork(8).Data;
    sum_sq_e_a    = block.Dwork(9).Data;
    mean_y_a      = block.Dwork(10).Data;
    M2_y_a        = block.Dwork(11).Data;
    nis_in        = block.Dwork(12).Data;
    nis_total     = block.Dwork(13).Data;

    % Trajectory metrics
    traj_SSE      = block.Dwork(14).Data;
    traj_SAE      = block.Dwork(15).Data;
    traj_N        = block.Dwork(16).Data;
    traj_sumY     = block.Dwork(17).Data;
    traj_sumY2    = block.Dwork(18).Data;

    % ---------- parameters ----------
    d = 0.9;
    p.A_xsec = pi/4*d^2;
    p.gamma0 = 0;
    p.gamma1 = 1;
    p.Y_O2_LAC_growth = 0.52;
    p.Y_O2_LAC_prod   = 0.12;

    % ---------- known feed ----------
    u = [88; 16; 3; 0.4; 60; 60];

    % ---------- process noise ----------
    % Treat these as continuous-time noise intensities per hour.
    Qc = diag([1e-2 1e-3 1e-3 1e-4 1e-4 1e-6]);
    Qd = Qc * dt;

    % ---------- measurement noise ----------
    sigma_V   = 10;
    sigma_CX  = 0.5;
    sigma_DO  = 5e-4;
    sigma_MEV = 0.01;   % g/L

    R_V   = sigma_V^2;
    R_CX  = sigma_CX^2;
    R_DO  = sigma_DO^2;
    R_MEV = sigma_MEV^2;

    % ---------- measurement models ----------
    H_V   = [1 0 0 0 0 0];
    H_CX  = [0 1 0 0 0 0];
    H_DO  = [0 0 0 0 0 1];
    H_MEV = [0 0 0 0 1 0];

    % ---------- measurement cadences in HOURS ----------
    Ts_V   = 5  / 3600;
    Ts_CX  = 10 / 3600;
    Ts_DO  = 5  / 3600;
    Ts_MEV = 24;

    % ---------- prediction ----------
    f = CSTRO_ode(xhat, u, p);
    xhat = xhat + dt * f;

    [A,~,~,~,~,~] = linearize_CSTRO(xhat, u, p);
    Ad = expm(A * dt);

    P = Ad * P * Ad' + Qd;
    P = 0.5 * (P + P');   % enforce symmetry

    % Save prior for diagnostics / NIS
    xhat_pred = xhat;
    P_pred    = P;

    % ---------- multi-rate update using ACTUAL TIME ----------
    Hk = [];
    yk = [];
    rdiag = [];

    mev_included = false;
    idx_mev = NaN;

    [take_V, next_V_time] = sample_due(t, next_V_time, Ts_V);
    if take_V
        Hk = [Hk; H_V];
        yk = [yk; y(1)];
        rdiag = [rdiag; R_V];
    end

    [take_CX, next_CX_time] = sample_due(t, next_CX_time, Ts_CX);
    if take_CX
        Hk = [Hk; H_CX];
        yk = [yk; y(2)];
        rdiag = [rdiag; R_CX];
    end

    [take_DO, next_DO_time] = sample_due(t, next_DO_time, Ts_DO);
    if take_DO
        Hk = [Hk; H_DO];
        yk = [yk; y(3)];
        rdiag = [rdiag; R_DO];
    end

    [take_MEV, next_MEV_time] = sample_due(t, next_MEV_time, Ts_MEV);
    if take_MEV
        idx_mev = size(Hk,1) + 1;
        Hk = [Hk; H_MEV];
        yk = [yk; y(4)];
        rdiag = [rdiag; R_MEV];
        mev_included = true;
    end

    K = zeros(6, max(1,size(Hk,1)));  % for safe indexing in diagnostics

    if ~isempty(Hk)
        Rk = diag(rdiag);

        innov = yk - Hk * xhat;
        S = Hk * P * Hk' + Rk;
        S = 0.5 * (S + S');

        K = (P * Hk') / S;

        xhat = xhat + K * innov;

        % Joseph stabilized covariance update
        I = eye(6);
        P = (I - K*Hk) * P * (I - K*Hk)' + K * Rk * K';
        P = 0.5 * (P + P');
    end

    % ---------- assay metrics / diagnostics ----------
    if mev_included
        innov_mev = y(4) - (H_MEV * xhat_pred);
        S_mev     = H_MEV * P_pred * H_MEV' + R_MEV;
        nis       = (innov_mev * innov_mev) / max(S_mev, eps);

        % 95% interval for chi-square(dof=1)
        nis_lo = 0.000982;
        nis_hi = 5.023886;

        nis_total = nis_total + 1;
        if (nis >= nis_lo) && (nis <= nis_hi)
            nis_in = nis_in + 1;
        end

        % assay error after update
        e_assay = y(4) - xhat(5);

        n_assay     = n_assay + 1;
        sum_abs_e_a = sum_abs_e_a + abs(e_assay);
        sum_sq_e_a  = sum_sq_e_a  + e_assay^2;

        % Welford update for assay variance
        if n_assay == 1
            mean_y_a = y(4);
            M2_y_a   = 0;
        else
            dy       = y(4) - mean_y_a;
            mean_y_a = mean_y_a + dy / n_assay;
            dy2      = y(4) - mean_y_a;
            M2_y_a   = M2_y_a + dy * dy2;
        end

        rmse_assay = sqrt(sum_sq_e_a / max(n_assay,1));
        mae_assay  = sum_abs_e_a / max(n_assay,1);
        sst_assay  = M2_y_a;
        if sst_assay > 1e-12
            R2_assay = 1 - (sum_sq_e_a / sst_assay);
        else
            R2_assay = NaN;
        end
        nis_cov = nis_in / max(nis_total,1);

        assignin('base', 'L1_assay_metrics', [rmse_assay, mae_assay, R2_assay, nis_cov]);

        % Log one row per assay:
        % [t_h, y_mev, xpred_mev, xpost_mev, innov_mev, P55_pred, K55_mev, nis]
        K55_mev = K(5, idx_mev);
        assay_row = [t, y(4), xhat_pred(5), xhat(5), innov_mev, P_pred(5,5), K55_mev, nis];

        try
            logM = evalin('base', 'L1_assay_log;');
            logM = [logM; assay_row];
            assignin('base', 'L1_assay_log', logM);
        catch
            assignin('base', 'L1_assay_log', assay_row);
        end
    end

    % ---------- full-trajectory metrics against truth ----------
    lov_true = x_true(5);
    lov_hat  = xhat(5);
    e_traj   = lov_true - lov_hat;

    traj_SSE   = traj_SSE + e_traj^2;
    traj_SAE   = traj_SAE + abs(e_traj);
    traj_N     = traj_N + 1;
    traj_sumY  = traj_sumY  + lov_true;
    traj_sumY2 = traj_sumY2 + lov_true^2;

    rmse_traj = sqrt(traj_SSE / max(traj_N,1));
    mae_traj  = traj_SAE / max(traj_N,1);

    meanY_traj = traj_sumY / max(traj_N,1);
    SST_traj   = traj_sumY2 - traj_N * meanY_traj^2;

    if SST_traj > 1e-12
        R2_traj = 1 - (traj_SSE / SST_traj);
    else
        R2_traj = NaN;
    end

    assignin('base', 'L1_traj_metrics', [rmse_traj, mae_traj, R2_traj]);

    % ---------- save back ----------
    block.Dwork(1).Data  = xhat;
    block.Dwork(2).Data  = reshape(P, 36, 1);

    block.Dwork(3).Data  = next_V_time;
    block.Dwork(4).Data  = next_CX_time;
    block.Dwork(5).Data  = next_DO_time;
    block.Dwork(6).Data  = next_MEV_time;

    block.Dwork(7).Data  = n_assay;
    block.Dwork(8).Data  = sum_abs_e_a;
    block.Dwork(9).Data  = sum_sq_e_a;
    block.Dwork(10).Data = mean_y_a;
    block.Dwork(11).Data = M2_y_a;
    block.Dwork(12).Data = nis_in;
    block.Dwork(13).Data = nis_total;

    block.Dwork(14).Data = traj_SSE;
    block.Dwork(15).Data = traj_SAE;
    block.Dwork(16).Data = traj_N;
    block.Dwork(17).Data = traj_sumY;
    block.Dwork(18).Data = traj_sumY2;
end

function [due, next_t] = sample_due(t, next_t, Ts)
% Returns true if a sample is due at current time t.
% Advances next_t to the first future sample time after t.
    tol = 1e-10;
    due = (t + tol >= next_t);

    if due
        nskip = floor((t - next_t) / Ts) + 1;
        if nskip < 1
            nskip = 1;
        end
        next_t = next_t + nskip * Ts;
    end
end