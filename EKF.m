function [xhat_out, yhat_out] = EKF_CSTR_baseline(u, y_meas, dt)
%#codegen
% EKF for KTB1 CSTR with baseline sensors: [h, C_DO, biomass-proxy]
% Inputs:
%   u      (6x1) : [C_X_F, C_LAC_F, C_N_F, C_MEV_F, Q_IN, Q_OUT]
%   y_meas (3x1) : [h, C_DO, C_X_proxy] from the plant (CSTRO)
%   dt     (1x1) : sample time [s]
% Outputs:
%   xhat_out (6x1) : EKF state estimate [V, C_X, C_LAC, C_N, C_MEV, C_DO]
%   yhat_out (3x1) : predicted measurements [h, C_DO, C_X_proxy]

% ---------- persistent filter state ----------
persistent p xhat P Q R gamma0 gamma1
if isempty(xhat)
    % params
    p = CSTR_default_params();
    % choose measurement set & proxy mapping
    p.meas_outputs = {'h','C_DO','C_X_proxy'};
    gamma0 = p.gamma0; gamma1 = p.gamma1;  % (set both in params if you like)

    % initial state guess (close to CSTRO IC; adjust if needed)
    xhat = [5000; 106.35; 14.81; 2.52; 1.20; 6e-3];

    % initial covariance
    P = diag([1e2, 1, 1, 1, 1, 1])*1e-2;

    % process/measurement noise (tunable; units: per hour model, seconds for dt)
    Q = diag([  (1e-3)^2, (5e-3)^2, (1e-2)^2, (5e-3)^2, (1e-3)^2, (5e-5)^2 ]);
    % baseline sensors: [h, C_DO, C_X_proxy]
    R = diag([ (2e-3)^2, (2e-4)^2, (5e-3)^2 ]);
end

% convert dt to hours for the continuous-time model
dt_hr = dt/3600;

% ---------- PREDICT ----------
[A,~,~,~,~,~] = linearize_CSTR(xhat, u, p);
f = CSTR_f(xhat, u, p);

% state prediction (Euler); replace with ode1/2 if you prefer
xhat = xhat + dt_hr * f;

% covariance prediction: Ad ~ I + dt*A ; Qd ~ Q*dt
Ad = eye(6) + dt_hr*A;
Qd = Q * dt_hr;
P  = Ad*P*Ad' + Qd;

% ---------- MEASUREMENT PREDICTION ----------
% y = [h, C_DO, C_X_proxy]
yhat = zeros(3,1);
V=xhat(1); C_X=xhat(2); C_DO=xhat(6);
h = (V/1000)/p.A_xsec;
yhat(1) = h;
yhat(2) = C_DO;
yhat(3) = gamma0 + gamma1*C_X;

% Measurement Jacobian H = dh/dx (3x6) via numeric Jacobian from linearize_CSTR
% (We can reuse C from linearize_CSTR if p.meas_outputs matches)
p.meas_outputs = {'h','C_DO','C_X_proxy'};
[~,~,H,~,~,~] = linearize_CSTR(xhat, u, p);

% ---------- UPDATE ----------
nu = y_meas - yhat;             % innovation
S  = H*P*H' + R;                % innovation covariance
K  = P*H' / S;                  % Kalman gain

xhat = xhat + K*nu;

% Joseph form for numerical stability
I = eye(6);
P = (I - K*H)*P*(I - K*H)' + K*R*K';

% outputs
xhat_out = xhat;
yhat_out = yhat;
end
