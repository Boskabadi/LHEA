function f = CSTRO_ode(x,u,p)
% States:
% x(1) V        [L]
% x(2) C_X      [g/L]
% x(3) C_LAC    [g/L]
% x(4) C_N      [g/L]
% x(5) C_MEV    [g/L]
% x(6) C_DO     [g/L]

% Inputs:
% u(1) C_X_F    [g/L]
% u(2) C_LAC_F  [g/L]
% u(3) C_N_F    [g/L]
% u(4) C_MEV_F  [g/L]
% u(5) Q_IN     [L/h]
% u(6) Q_OUT    [L/h]

%% States
V        = x(1);
C_X      = x(2);
C_LAC    = x(3);
C_N      = x(4);
C_MEV    = x(5);
C_DO     = x(6);

%% Inputs
C_X_F    = u(1);
C_LAC_F  = u(2);
C_N_F    = u(3);
C_MEV_F  = u(4);
Q_IN     = u(5);
Q_OUT    = u(6);

%% Kinetic parameters (units preserved)
mu_max     = 0.120;      % [1/h]
q_max_MEV  = 1.9e-4;     % [g_MEV/(g_X·h)]
K_MEV      = 5.42e-6;    % [-]
Y_X_LAC    = 0.483;      % [g_X/g_LAC]
Y_X_N      = 20.06;      % [g_X/g_N]
Y_MEV_LAC  = 7.06e-4;    % [g_MEV/g_LAC]

K_LAC      = 1.63;       % [g_LAC/g_X]
K_N        = 8.84e-2;    % [g_N/g_X]
K_LAC_MEV  = 13.23;      % [g_LAC/g_X]
K_I_N      = 0.158;      % [g_N/L]
K_I_N_MEV  = 9.65e-2;    % [g_N/L]

%% Specific rates
mu_1 = mu_max * (C_LAC/(C_LAC + K_LAC*C_X)) ...
               * (C_N  /(C_N  + K_N*C_X));

mu_2 = mu_1 * (K_I_N/(K_I_N + C_N));

mu_3 = q_max_MEV * (C_LAC/(C_LAC + K_LAC_MEV*C_X)) ...
                   * (K_I_N_MEV/(K_I_N_MEV + C_N));

%% Yield conversions
Y_LAC_X   = 1/Y_X_LAC;      % [g_LAC/g_X]
Y_LAC_MEV = 1/Y_MEV_LAC;    % [g_LAC/g_MEV]
Y_N_X     = 1/Y_X_N;        % [g_N/g_X]

%% Oxygen parameters
Kla      = 120;       % [1/h]
C_DO_sat = 0.008;     % [g/L]
C_DO_F   = 0.0;       % [g/L]

% OUR coefficients (lumped, unit-consistent)
Y_O2_LAC_growth = p.Y_O2_LAC_growth; % [g_O2/g_LAC]
Y_O2_LAC_prod   = p.Y_O2_LAC_prod;   % [g_O2/g_LAC]

%% Volume balance
dVdt = Q_IN - Q_OUT;     % [L/h]

%% Biomass
dC_Xdt = mu_1*C_X ...
       + (Q_IN/V)*(C_X_F - C_X);

%% Lactose
dC_LACdt = -(Y_LAC_X*mu_2 + Y_LAC_MEV*mu_3)*C_X ...
         + (Q_IN/V)*(C_LAC_F - C_LAC);

%% Nitrogen
dC_Ndt = -(Y_N_X*mu_1)*C_X ...
       + (Q_IN/V)*(C_N_F - C_N);

%% Product
dC_MEVdt = (mu_3 + K_MEV*C_LAC)*C_X ...
         + (Q_IN/V)*(C_MEV_F - C_MEV);

%% Oxygen
q_lac_growth = Y_LAC_X   * mu_2;     % [g_LAC/g_X/h]
q_lac_prod   = Y_LAC_MEV * mu_3;     % [g_LAC/g_X/h]

OUR = (q_lac_growth*Y_O2_LAC_growth ...
     + q_lac_prod*Y_O2_LAC_prod) * C_X;   % [g_O2/L/h]

OTR = Kla*(C_DO_sat - C_DO);               % [g_O2/L/h]

dC_DOdt = OTR - OUR + (Q_IN/V)*(C_DO_F - C_DO);

%% Assemble
f = [dVdt;
     dC_Xdt;
     dC_LACdt;
     dC_Ndt;
     dC_MEVdt;
     dC_DOdt];
end
