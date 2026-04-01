function EKF_CSTR_L33(block)
% L3 further-improved hybrid EKF
%
% Improvements over L32:
% 1) split MEV correction into:
%    - slow bias rate d_bias
%    - fast feature-based residual theta' * phi
% 2) learn on assay-interval average features/rates for better conditioning
% 3) allow slightly stronger MEV model mismatch through Q and correction bounds
%
% Inputs:
%   1) x_true (6)   : only for simulation metrics
%   2) y_meas (4)   : [V; CX; DO; MEV_assay]
%   3) dt (1) [h]
%
% Output:
%   [xhat(6); qscale]  -> 7x1
%
% IMPORTANT:
% This code assumes block.CurrentTime is in HOURS.

    setup(block);
end

function setup(block)
    block.NumInputPorts  = 3;
    block.NumOutputPorts = 1;

    block.InputPort(1).Dimensions = 6;   % x_true
    block.InputPort(2).Dimensions = 4;   % [V; CX; DO; MEV]
    block.InputPort(3).Dimensions = 1;   % dt [h]

    for i = 1:3
        block.InputPort(i).DatatypeID = 0;
        block.InputPort(i).DirectFeedthrough = true;
    end

    block.OutputPort(1).Dimensions = 7;  % 6 states + qscale
    block.OutputPort(1).DatatypeID = 0;

    block.SampleTimes = [0 0];
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('PostPropagationSetup', @PostPropSetup);
    block.RegBlockMethod('InitializeConditions', @InitConditions);
    block.RegBlockMethod('Outputs',              @Output);
    block.RegBlockMethod('Update',               @Update);
end

function PostPropSetup(block)
    % Dworks:
    %  1  xhat(6)
    %  2  P(36)
    %  3  next_V_time
    %  4  next_CX_time
    %  5  next_DO_time
    %  6  next_MEV_time
    %  7  qscale
    %  8  theta(4)
    %  9  Ptheta(16)
    % 10  phi_int(4)
    % 11  d_bias
    % 12  t_acc
    %
    % assay metrics
    % 13  n_assay
    % 14  sum_abs_err_assay
    % 15  sum_sq_err_assay
    % 16  mean_y_assay
    % 17  M2_y_assay
    % 18  nis_in
    % 19  nis_total
    %
    % full trajectory metrics
    % 20  traj_SSE
    % 21  traj_SAE
    % 22  traj_N
    % 23  traj_sumY
    % 24  traj_sumY2

    block.NumDworks = 24;

    for i = 1:24
        block.Dwork(i).DatatypeID = 0;
        block.Dwork(i).Complexity = 'Real';
        block.Dwork(i).UsedAsDiscState = true;
    end

    block.Dwork(1).Name  = 'xhat';               block.Dwork(1).Dimensions  = 6;
    block.Dwork(2).Name  = 'P';                  block.Dwork(2).Dimensions  = 36;
    block.Dwork(3).Name  = 'next_V_time';        block.Dwork(3).Dimensions  = 1;
    block.Dwork(4).Name  = 'next_CX_time';       block.Dwork(4).Dimensions  = 1;
    block.Dwork(5).Name  = 'next_DO_time';       block.Dwork(5).Dimensions  = 1;
    block.Dwork(6).Name  = 'next_MEV_time';      block.Dwork(6).Dimensions  = 1;
    block.Dwork(7).Name  = 'qscale';             block.Dwork(7).Dimensions  = 1;
    block.Dwork(8).Name  = 'theta';              block.Dwork(8).Dimensions  = 4;
    block.Dwork(9).Name  = 'Ptheta';             block.Dwork(9).Dimensions  = 16;
    block.Dwork(10).Name = 'phi_int';            block.Dwork(10).Dimensions = 4;
    block.Dwork(11).Name = 'd_bias';             block.Dwork(11).Dimensions = 1;
    block.Dwork(12).Name = 't_acc';              block.Dwork(12).Dimensions = 1;

    block.Dwork(13).Name = 'n_assay';            block.Dwork(13).Dimensions = 1;
    block.Dwork(14).Name = 'sum_abs_err_assay';  block.Dwork(14).Dimensions = 1;
    block.Dwork(15).Name = 'sum_sq_err_assay';   block.Dwork(15).Dimensions = 1;
    block.Dwork(16).Name = 'mean_y_assay';       block.Dwork(16).Dimensions = 1;
    block.Dwork(17).Name = 'M2_y_assay';         block.Dwork(17).Dimensions = 1;
    block.Dwork(18).Name = 'nis_in';             block.Dwork(18).Dimensions = 1;
    block.Dwork(19).Name = 'nis_total';          block.Dwork(19).Dimensions = 1;

    block.Dwork(20).Name = 'traj_SSE';           block.Dwork(20).Dimensions = 1;
    block.Dwork(21).Name = 'traj_SAE';           block.Dwork(21).Dimensions = 1;
    block.Dwork(22).Name = 'traj_N';             block.Dwork(22).Dimensions = 1;
    block.Dwork(23).Name = 'traj_sumY';          block.Dwork(23).Dimensions = 1;
    block.Dwork(24).Name = 'traj_sumY2';         block.Dwork(24).Dimensions = 1;
end

function InitConditions(block)
    x0 = [5000; 106; 14.8; 2.52; 1.2; 0.008];
    P0 = diag([1e2 5 5 5 5 1e-4]);

    block.Dwork(1).Data = x0;
    block.Dwork(2).Data = reshape(P0,36,1);

    block.Dwork(3).Data = 0;
    block.Dwork(4).Data = 0;
    block.Dwork(5).Data = 0;
    block.Dwork(6).Data = 0;

    block.Dwork(7).Data  = 1.0;                    % qscale
    block.Dwork(8).Data  = zeros(4,1);            % theta
    block.Dwork(9).Data  = reshape(25*eye(4),16,1); % Ptheta
    block.Dwork(10).Data = zeros(4,1);            % phi_int
    block.Dwork(11).Data = 0.0;                   % d_bias [g/L/h]
    block.Dwork(12).Data = 0.0;                   % accumulated time since last assay [h]

    block.Dwork(13).Data = 0;
    block.Dwork(14).Data = 0;
    block.Dwork(15).Data = 0;
    block.Dwork(16).Data = 0;
    block.Dwork(17).Data = 0;
    block.Dwork(18).Data = 0;
    block.Dwork(19).Data = 0;

    block.Dwork(20).Data = 0;
    block.Dwork(21).Data = 0;
    block.Dwork(22).Data = 0;
    block.Dwork(23).Data = 0;
    block.Dwork(24).Data = 0;

    assignin('base','L3_assay_log',zeros(0,17));
    assignin('base','L3_assay_log_cols', ...
        {'t_h','y_mev','xpred_mev','xpost_mev','innov_mev','P55_pred', ...
         'K55_mev','nis','qscale','target_rate','d_bias','r_feat_pred', ...
         'r_total_pred','theta_1','theta_2','theta_3','theta_4'});
    assignin('base','L3_assay_metrics',[NaN NaN NaN NaN]);
    assignin('base','L3_traj_metrics',[NaN NaN NaN]);
    assignin('base','L3_metrics',[NaN NaN NaN NaN]);
    assignin('base','L3_theta',zeros(4,1));
    assignin('base','L3_qscale',1.0);
    assignin('base','L3_dbias',0.0);
end

function Output(block)
    xhat   = block.Dwork(1).Data;
    qscale = block.Dwork(7).Data;
    block.OutputPort(1).Data = [xhat; qscale];
end

function Update(block)
    x_true = block.InputPort(1).Data;
    y      = block.InputPort(2).Data;    % [V; CX; DO; MEV_assay]
    dt     = block.InputPort(3).Data;
    dt     = max(dt, eps);

    t = block.CurrentTime;

    xhat = block.Dwork(1).Data;
    P    = reshape(block.Dwork(2).Data,6,6);

    next_V_time   = block.Dwork(3).Data;
    next_CX_time  = block.Dwork(4).Data;
    next_DO_time  = block.Dwork(5).Data;
    next_MEV_time = block.Dwork(6).Data;

    qscale  = block.Dwork(7).Data;
    theta   = block.Dwork(8).Data;
    Ptheta  = reshape(block.Dwork(9).Data,4,4);
    phi_int = block.Dwork(10).Data;
    d_bias  = block.Dwork(11).Data;
    t_acc   = block.Dwork(12).Data;

    n_assay     = block.Dwork(13).Data;
    sum_abs_e_a = block.Dwork(14).Data;
    sum_sq_e_a  = block.Dwork(15).Data;
    mean_y_a    = block.Dwork(16).Data;
    M2_y_a      = block.Dwork(17).Data;
    nis_in      = block.Dwork(18).Data;
    nis_total   = block.Dwork(19).Data;

    traj_SSE    = block.Dwork(20).Data;
    traj_SAE    = block.Dwork(21).Data;
    traj_N      = block.Dwork(22).Data;
    traj_sumY   = block.Dwork(23).Data;
    traj_sumY2  = block.Dwork(24).Data;

    % ---------- parameters ----------
    d = 0.9;
    p.A_xsec = pi/4*d^2;
    p.gamma0 = 0;
    p.gamma1 = 1;
    p.Y_O2_LAC_growth = 0.52;
    p.Y_O2_LAC_prod   = 0.12;
    p.q_max_scale     = qscale;

    % ---------- known feed ----------
    u = [88; 16; 3.0; 0.4; 60; 60];

    % ---------- process noise ----------
    % More freedom on MEV state than L31/L32
    Qc = diag([1e-2 1e-3 1e-3 1e-4 1.5e-3 1e-6]);
    Qd = Qc * dt;

    % ---------- measurement noise ----------
    sigma_V   = 10;
    sigma_CX  = 0.5;
    sigma_DO  = 5e-4;
    sigma_MEV = 0.01;   % slightly stronger trust in assay

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
    Ts_MEV = 24;   % keep identical to experiment and to L1 comparison

    % ---------- prediction ----------
    phi = residual_features(xhat);
    phi_int = phi_int + dt * phi;
    t_acc   = t_acc + dt;

    d_bias_max  = 2.5e-3;   % slow bias rate bound [g/L/h]
    r_feat_max  = 2.5e-3;   % fast feature residual bound [g/L/h]
    r_total_max = 4.5e-3;   % total correction bound [g/L/h]

    d_bias = clip_scalar(d_bias, -d_bias_max, d_bias_max);

    r_feat_pred = theta.' * phi;
    r_feat_pred = clip_scalar(r_feat_pred, -r_feat_max, r_feat_max);

    r_total_pred = d_bias + r_feat_pred;
    r_total_pred = clip_scalar(r_total_pred, -r_total_max, r_total_max);

    f = CSTRO_ode(xhat, u, p);
    f(5) = f(5) + r_total_pred;

    xhat = xhat + dt * f;
    xhat = clamp_state(xhat);

    [A,~,~,~,~,~] = linearize_CSTRO(xhat, u, p);
    Ad = expm(A * dt);

    P = Ad * P * Ad' + Qd;
    P = 0.5 * (P + P');

    xhat_pred = xhat;
    P_pred    = P;

    % ---------- multi-rate update ----------
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

    K = zeros(6, max(1,size(Hk,1)));

    if ~isempty(Hk)
        Rk = diag(rdiag);

        innov = yk - Hk * xhat;
        S = Hk * P * Hk' + Rk;
        S = 0.5 * (S + S');

        K = (P * Hk') / S;

        xhat = xhat + K * innov;
        xhat = clamp_state(xhat);

        I = eye(6);
        P = (I - K*Hk) * P * (I - K*Hk)' + K * Rk * K';
        P = 0.5 * (P + P');
    end

    % ---------- assay-time learning ----------
    if mev_included
        innov_mev = y(4) - (H_MEV * xhat_pred);
        S_mev     = H_MEV * P_pred * H_MEV' + R_MEV;
        nis       = (innov_mev * innov_mev) / max(S_mev, eps);

        nis_lo = 0.000982;
        nis_hi = 5.023886;

        nis_total = nis_total + 1;
        if (nis >= nis_lo) && (nis <= nis_hi)
            nis_in = nis_in + 1;
        end

        t_eff = max(t_acc, dt);
        zbar  = phi_int / t_eff;          % average feature over assay interval
        zbar  = max(min(zbar, 3), -3);

        target_rate = innov_mev / t_eff;  % rate error [g/L/h]

        % learning gate
        learn_ok = isfinite(nis) && (nis < 16) && all(isfinite(zbar));

        if learn_ok
            % ---- slow bias learner ----
            pred_rate = d_bias + theta.' * zbar;
            err_rate  = target_rate - pred_rate;

            eta_b = 0.18;   % slow but meaningful
            d_bias = d_bias + eta_b * err_rate;
            d_bias = clip_scalar(d_bias, -d_bias_max, d_bias_max);

            % ---- fast feature learner (RLS) ----
            % learn residual after bias update
            pred_rate2 = d_bias + theta.' * zbar;
            err_rate2  = target_rate - pred_rate2;

            lambda = 0.992;
            denom  = lambda + (zbar.' * Ptheta * zbar);
            Kth    = (Ptheta * zbar) / denom;

            theta  = theta + Kth * err_rate2;
            Ptheta = (Ptheta - Kth * (zbar.' * Ptheta)) / lambda;
            Ptheta = 0.5 * (Ptheta + Ptheta.');

            % enforce instantaneous residual limit
            phi_now = residual_features(xhat);
            r_feat_now = theta.' * phi_now;
            if abs(r_feat_now) > r_feat_max
                theta = theta * (r_feat_max / abs(r_feat_now));
            end

            % very slow qscale adaptation from the low-frequency bias only
            eta_q = 0.02;
            qscale = qscale + eta_q * d_bias;
            qscale = min(max(qscale, 0.7), 1.3);
        end

        % reset interval accumulators after assay
        phi_int = zeros(4,1);
        t_acc   = 0.0;

        % ---------- assay metrics ----------
        e_assay = y(4) - xhat(5);

        n_assay     = n_assay + 1;
        sum_abs_e_a = sum_abs_e_a + abs(e_assay);
        sum_sq_e_a  = sum_sq_e_a  + e_assay^2;

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

        assignin('base','L3_assay_metrics',[rmse_assay, mae_assay, R2_assay, nis_cov]);
        assignin('base','L3_metrics',[rmse_assay, mae_assay, R2_assay, nis_cov]);
        assignin('base','L3_theta',theta);
        assignin('base','L3_qscale',qscale);
        assignin('base','L3_dbias',d_bias);

        K55_mev = K(5, idx_mev);

        assay_row = [t, y(4), xhat_pred(5), xhat(5), innov_mev, P_pred(5,5), ...
                     K55_mev, nis, qscale, target_rate, d_bias, r_feat_pred, ...
                     r_total_pred, theta(1), theta(2), theta(3), theta(4)];

        try
            logM = evalin('base','L3_assay_log;');
            logM = [logM; assay_row];
            assignin('base','L3_assay_log',logM);
        catch
            assignin('base','L3_assay_log',assay_row);
        end
    end

    % ---------- full-trajectory metrics ----------
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

    assignin('base','L3_traj_metrics',[rmse_traj, mae_traj, R2_traj]);

    % ---------- save back ----------
    block.Dwork(1).Data  = xhat;
    block.Dwork(2).Data  = reshape(P,36,1);

    block.Dwork(3).Data  = next_V_time;
    block.Dwork(4).Data  = next_CX_time;
    block.Dwork(5).Data  = next_DO_time;
    block.Dwork(6).Data  = next_MEV_time;

    block.Dwork(7).Data  = qscale;
    block.Dwork(8).Data  = theta;
    block.Dwork(9).Data  = reshape(Ptheta,16,1);
    block.Dwork(10).Data = phi_int;
    block.Dwork(11).Data = d_bias;
    block.Dwork(12).Data = t_acc;

    block.Dwork(13).Data = n_assay;
    block.Dwork(14).Data = sum_abs_e_a;
    block.Dwork(15).Data = sum_sq_e_a;
    block.Dwork(16).Data = mean_y_a;
    block.Dwork(17).Data = M2_y_a;
    block.Dwork(18).Data = nis_in;
    block.Dwork(19).Data = nis_total;

    block.Dwork(20).Data = traj_SSE;
    block.Dwork(21).Data = traj_SAE;
    block.Dwork(22).Data = traj_N;
    block.Dwork(23).Data = traj_sumY;
    block.Dwork(24).Data = traj_sumY2;
end

function [due, next_t] = sample_due(t, next_t, Ts)
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

function phi = residual_features(x)
% Operating-point features for fast residual learner
% x = [V; CX; LAC; N; MEV; DO]
    zCX  = (x(2) - 100.0) / 20.0;
    zN   = (x(4) - 2.50)  / 0.40;
    zDO  = (x(6) - 0.008) / 0.002;
    zMEV = (x(5) - 1.20)  / 0.12;

    phi = [zCX; zN; zDO; zMEV];
    phi = max(min(phi, 2.5), -2.5);
end

function x = clamp_state(x)
    x(1) = max(x(1), 1.0);      % V
    x(2) = max(x(2), 0.0);      % CX
    x(3) = max(x(3), 0.0);      % LAC
    x(4) = max(x(4), 0.0);      % N
    x(5) = max(x(5), 0.0);      % MEV
    x(6) = max(x(6), 1e-8);     % DO
end

function y = clip_scalar(x, lo, hi)
    y = min(max(x, lo), hi);
end