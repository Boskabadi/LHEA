function p = CSTR_default_params()
% Codegen-friendly parameter struct for KTB1 CSTR (no cell fields).

p = struct( ...
  'mu_max',0,'q_max_MEV',0,'K_MEV',0, ...
  'Y_X_LAC',0,'Y_X_N',0,'Y_MEV_LAC',0, ...
  'K_LAC',0,'K_N',0,'K_LAC_MEV',0,'K_I_N',0,'K_I_N_MEV',0, ...
  'd_reactor',0,'A_xsec',0,'rhoL',1.0, ...
  'Kla',0,'CDO_sat',0,'CDO_F',0, ...
  'Y_O2_Lac_growth',0,'Y_O2_Lac_prod',0, ...
  'gamma0',0,'gamma1',1, ...
  'eps',1e-12);

% ---- Kinetics ----
p.mu_max    = 0.120;      % 1/h
p.q_max_MEV = 1.9e-4;     % g_MEV/g_X/h
p.K_MEV     = 5.42e-6;
p.Y_X_LAC   = 0.483;      % g_X/g_LAC
p.Y_X_N     = 20.06;      % g_X/g_N
p.Y_MEV_LAC = 7.06e-4;    % g_MEV/g_LAC
p.K_LAC     = 1.63;
p.K_N       = 8.84e-2;
p.K_LAC_MEV = 13.23;
p.K_I_N     = 0.158;
p.K_I_N_MEV = 9.65e-2;

% ---- Geometry ----
p.d_reactor = 0.9;                       % m
p.A_xsec    = pi/4 * p.d_reactor^2;      % m^2

% ---- Oxygen ----
p.Kla     = 120;         % 1/h
p.CDO_sat = 8e-3;        % g/L
p.CDO_F   = 0.0;         % g/L

% ---- O2 yields (your stoichiometry) ----
Yxs   = p.Y_X_LAC;  Yps = p.Y_MEV_LAC;
MW_S  = 342.297;  OX=0.52; HX=1.81; NX=0.21;
MW_XC = 12 + HX + 16*OX + 14*NX;  MW_P = 404.547;  MW_O2 = 32;

b = Yxs*MW_S/(MW_XC*12);
a = (OX*b)/2 + 1 - (HX*b)/4 - b/4;
p.Y_O2_Lac_growth = a*MW_O2/MW_S;

f = Yps*MW_S/MW_P;  g = 12 - 24*f;  h2o = (22 - 36*f)/2;
e = ((5*f + 2*g + h2o) - 11)/2;
p.Y_O2_Lac_prod = e*MW_O2/MW_S;

% ---- Biomass proxy mapping ----
p.gamma0 = 0;  % offset
p.gamma1 = 1;  % slope (proxy = gamma0 + gamma1*C_X)
end
