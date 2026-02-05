function [sys,X_0,str,ts] = CSTRO(t,x,u,flag,X_0)

switch flag
    case 0 % initialization
        sizes = simsizes;
        sizes.NumContStates  = 6;   % V, C_X, C_LAC, C_N, C_MEV, C_DO
        sizes.NumDiscStates  = 0;
        sizes.NumOutputs     = 7;   % V, C_X, C_LAC, C_N, C_MEV, C_DO, h
        sizes.NumInputs      = 6;   % C_X_F, C_LAC_F, C_N_F, C_MEV_F, Q_IN, Q_OUT
        sizes.DirFeedthrough = 0;
        sizes.NumSampleTimes = 1;
        sys = simsizes(sizes);

        % Initial conditions (adjust as needed)
        % V [L], C_X [g/L], C_LAC [g/L], C_N [g/L], C_MEV [g/L], C_DO [g/L]
        X_0 = [5000; 106.351163614666; 14.8096563646109; 2.51553709810012; 1.20000000000003; 0.008]; 
        str = [];
        ts  = [0 0];

    case 1 % derivatives
        % ---- States ----
        V        = x(1);  % L
        C_X_S    = x(2);  % g_X/L
        C_LAC_S  = x(3);  % g_LAC/L
        C_N_S    = x(4);  % g_N/L (adenine)
        C_MEV_S  = x(5);  % g_MEV/L
        C_DO     = x(6);  % g_O2/L

        % ---- Inputs ----
        C_X_F    = u(1);  % g_X/L
        C_LAC_F  = u(2);  % g_LAC/L
        C_N_F    = u(3);  % g_N/L (adenine)
        C_MEV_F  = u(4);  % g_MEV/L
        Q_IN     = u(5);  % L/h
        Q_OUT    = u(6);  % L/h

        % ---- Kinetic parameters (as in your original model) ----
        mu_max     = 0.120;     % 1/h

        % ===== L3 TEST 1: TRUE q_max_MEV step at t = 200 h =====
        q_max_MEV_nom = 1.9e-4;                 % g_MEV/g_X/h
        if t >= 50
            q_max_MEV = 1.2 * q_max_MEV_nom;    % +20% after 50 h (plant mismatch)
        else
            q_max_MEV = q_max_MEV_nom;
        end
        % =======================================================

        K_MEV      = 5.42e-6;
        Y_X_LAC    = 0.483;
        Y_X_N      = 20.06;
        Y_MEV_LAC  = 7.06e-4;
        K_LAC      = 1.63;
        K_N        = 8.84e-2;
        K_LAC_MEV  = 13.23;
        K_I_N      = 0.158;
        K_I_N_MEV  = 9.65e-2;


        % ---- Geometry ----
        d  = 0.9;                   % m, reactor diameter
        A  = pi/4*d^2;              % m^2
        h  = (V/1000)/A;            % m, liquid height (V [L] -> m^3)

        % ---- Specific rates ----
        % Growth terms (keep your original structure)
        mu_1 = mu_max * (C_LAC_S/(C_LAC_S + K_LAC*C_X_S)) * (C_N_S/(C_N_S + K_N*C_X_S));
        mu_2 = mu_max * (C_LAC_S/(C_LAC_S + K_LAC*C_X_S)) * (C_N_S/(C_N_S + K_N*C_X_S)) * (K_I_N/(K_I_N + C_N_S));
        mu_3 = q_max_MEV * (C_LAC_S/(C_LAC_S + K_LAC_MEV*C_X_S)) * (K_I_N_MEV/(K_I_N_MEV + C_N_S)); % g_MEV/g_X/h

        % Yield conversions used in mass balances (your original)
        Y_LAC_X   = 1/Y_X_LAC;      % g_LAC / g_X
        Y_LAC_MEV = 1/Y_MEV_LAC;    % g_LAC / g_MEV
        Y_N_X     = 1/Y_X_N;        % g_N / g_X

        % ---- Stoichiometry for DO/CO2 (your provided section, cleaned) ----
        % Constants
        Yxs   = 0.483;      % g_X per g_LAC (same as Y_X_LAC)
        Yps   = 7.06e-4;    % g_MEV per g_LAC (same as Y_MEV_LAC)
        MW_S  = 342.297;    % g/mol lactose
        OX    = 0.52; HX = 1.81; NX = 0.21;           % C-mol composition of biomass
        MW_XC = 12 + HX + OX*16 + NX*14;              % g per C-mol biomass
        MW_P  = 404.547;    % g/mol (mevalonate or target product; per your note)
        MW_O2 = 32;         % g/mol O2
        MW_CO2= 44;         % g/mol CO2

        % Pathway 1: lactose -> biomass (+ CO2 + H2O)
        b = Yxs*MW_S/(MW_XC*12);                           % mol biomass per C-mol lactose
        a = (OX*b)/2 + 1 - (HX*b)/4 - b/4;                 % mol O2 per C-mol lactose
        c = 1 - b;                                         % mol CO2 per C-mol lactose
        dH2O = 0.92 + (3*b)/2 - (HX*b)/2;                  % mol H2O per C-mol lactose
        Y_O2_Lac_growth = a*MW_O2/MW_S;                    % g O2 per g lactose
        Y_CO2_Lac_growth = c*MW_CO2/MW_S;                  % g CO2 per g lactose (not used here)

        % Pathway 2: lactose -> product (+ CO2 + H2O)
        f = Yps*MW_S/MW_P;                                 % mol product per mol lactose
        g = 12 - f*24;                                     % mol CO2 per C-mol lactose
        h2o = (22 - f*36)/2;                               % mol H2O per C-mol lactose
        e = ((5*f + 2*g + h2o) - 11)/2;                    % mol O2 per C-mol lactose
        Y_O2_Lac_prod  = e*MW_O2/MW_S;                     % g O2 per g lactose
        % Y_CO2_Lac_prod = g*MW_CO2/MW_S;                  % g CO2 per g lactose (not used here)

        % ---- Oxygen transfer parameters ----
        Kla      = 120;        % 1/h (example; set per your system)
        C_DO_sat = 0.008;      % g/L (e.g., ~8 mg/L at 25°C, 1 atm; adjust as needed)
        C_DO_F   = 0.0;        % g/L in feed (often negligible)

        % ---- Mass balances ----
        dVdt        = (Q_IN - Q_OUT);

        % Biomass
        dC_X_Sdt    = (mu_1*C_X_S + (Q_IN*C_X_F - Q_OUT*C_X_S - C_X_S*(Q_IN - Q_OUT))/V);

        % Lactose
        dC_LAC_Sdt  = (-(Y_LAC_X*mu_2 + Y_LAC_MEV*mu_3)*C_X_S + (Q_IN*C_LAC_F - Q_OUT*C_LAC_S - C_LAC_S*(Q_IN - Q_OUT))/V);

        % Nitrogen (adenine)
        dC_N_Sdt    = (-(Y_N_X*mu_1)*C_X_S + (Q_IN*C_N_F - Q_OUT*C_N_S - C_N_S*(Q_IN - Q_OUT))/V);

        % Product
        dC_MEV_Sdt  = ((mu_3 + K_MEV*C_LAC_S)*C_X_S + (Q_IN*C_MEV_F - Q_OUT*C_MEV_S - C_MEV_S*(Q_IN - Q_OUT))/V);

        % ---- Oxygen uptake & transfer ----
        % Lactose uptake split by pathway [g_lac/g_X/h]
        q_lac_growth = Y_LAC_X   * mu_2;
        q_lac_prod   = Y_LAC_MEV * mu_3;

        % OUR [g_O2/(L h)]
        OUR = (q_lac_growth * Y_O2_Lac_growth + q_lac_prod * Y_O2_Lac_prod) * C_X_S;

        % OTR [g_O2/(L h)]
        OTR = Kla * (C_DO_sat - C_DO);

        % DO state
        dC_DO_dt = OTR - OUR + (Q_IN*C_DO_F - Q_OUT*C_DO - C_DO*(Q_IN - Q_OUT))/V;

        sys = [dVdt; dC_X_Sdt; dC_LAC_Sdt; dC_N_Sdt; dC_MEV_Sdt; dC_DO_dt];

    case 3 % outputs
        d = 0.9;  A = pi/4*d^2;
        h = (x(1)/1000)/A;   % m

        % Expose C_DO as well
        sys = [x(1); x(2); x(3); x(4); x(5); x(6); h];

    case {2,4,9}
        sys = [];
    otherwise
        error(['unhandled flag = ', num2str(flag)]);
end
