%% run_multirate_observability.m
% Multi-rate observability analysis for CSTRO
% Uses: CSTRO_ode.m, linearize_CSTRO.m

clear; clc; close all;

%% 1. Nominal operating point (same as in CSTRO.m)
% States: [V, C_X, C_LAC, C_N, C_MEV, C_DO]
x0 = [5000; ...       % V [L]
      106.3512; ...   % C_X [g/L]
      14.8097; ...    % C_LAC [g/L]
      2.5155; ...     % C_N [g/L]
      1.2000; ...     % C_MEV [g/L]
      0.0060];        % C_DO [g/L]

% Inputs: [C_X_F, C_LAC_F, C_N_F, C_MEV_F, Q_IN, Q_OUT]
u0 = [0; 0; 0; 0; 0; 0];

%% 2. Parameters for linearization and O2 model
% Geometry (same as CSTRO.m)
d = 0.9;                     % m, reactor diameter
p.A_xsec = pi/4*d^2;         % m^2, cross-sectional area

% Biomass proxy parameters (can be calibrated)
p.gamma0 = 0;
p.gamma1 = 1;                % proxy ≈ C_X

% Oxygen stoichiometry (example values – adjust if you have precise ones)
p.Y_O2_LAC_growth = 0.52;    % g O2 / g lactose (growth)
p.Y_O2_LAC_prod   = 0.12;    % g O2 / g lactose (product)

%% 3. Linearize the nonlinear model at (x0, u0)
[A,B,C_full,D,f0,y0] = linearize_CSTRO(x0,u0,p);
nx = size(A,1);

%% 4. Base discrete-time step and horizon
% Time unit in the model is hours.
% We choose dt = 1 s = 1/3600 h as the base step for the multi-rate grid.
dt_sec        = 1;                 % base step in seconds
dt_h          = dt_sec/3600;       % base step in hours
T_horizon_hr  = 1;                 % analyze observability over 1 hour
N_h           = round(T_horizon_hr*3600/dt_sec);   % number of base steps

fprintf('Base step: %.2f s, Horizon: %.2f h, Steps: %d\n', ...
        dt_sec, T_horizon_hr, N_h);

%% 5. Define realistic sensor cadences (in seconds)
% You can adapt these numbers to your setup
cad.DO_fast_sec       = 5;      % DO probe (1–5 s)
cad.Level_sec         = 5;      % Level sensor (1–5 s)
cad.Biomass_fast_sec  = 10;     % Capacitance (5–15 s)
cad.Lactose_Raman_sec = 60;     % Raman/NIR (30–120 s)
cad.MEV_HPLC_sec      = 1800;   % HPLC (30 min) – slow

%% 6. Define multi-rate measurement scenarios
% State indices: [1 V, 2 C_X, 3 C_LAC, 4 C_N, 5 C_MEV, 6 C_DO]
%
% Each scenario: { name, measured_state_indices, sampling_times_sec }
% sampling_times_sec must be same length as measured_state_indices

scenarios = {
    % 1) Only fast DO
    {'Senario 1', ...
        [6], ...
        [cad.DO_fast_sec]}

    % 2) Fast DO + level
    {'Senario 2', ...
        [6 1], ...
        [cad.DO_fast_sec cad.Level_sec]}

    % 3) DO + biomass (both fast)
    {'Senario 3', ...
        [6 2], ...
        [cad.DO_fast_sec cad.Biomass_fast_sec]}

    % 4) DO + biomass + level (all fast)
    {'Senario 4', ...
        [6 2 1], ...
        [cad.DO_fast_sec cad.Biomass_fast_sec cad.Level_sec]}

    % 5) DO fast + Lactose soft-sensor (60 s)
    {'Senario 5', ...
        [6 3], ...
        [cad.DO_fast_sec cad.Lactose_Raman_sec]}

    % 6) DO + biomass + Lactose soft-sensor
    {'Senario 6', ...
        [6 2 3], ...
        [cad.DO_fast_sec cad.Biomass_fast_sec cad.Lactose_Raman_sec]}

    % 7) Full realistic online/soft-sensor set (no MEV, no N)
    {'Senario 7', ...
        [1 6 2 3], ...
        [cad.Level_sec cad.DO_fast_sec cad.Biomass_fast_sec cad.Lactose_Raman_sec]}

    % 8) Add slow MEV HPLC (30 min)
    {'Senario 8', ...
        [1 6 2 3 5], ...
        [cad.Level_sec cad.DO_fast_sec cad.Biomass_fast_sec cad.Lactose_Raman_sec cad.MEV_HPLC_sec]}
};

nScen = numel(scenarios);
obs_ranks = zeros(nScen,1);

%% 7. Compute observability rank for each multi-rate scenario
fprintf('\n=== MULTI-RATE OBSERVABILITY ANALYSIS ===\n');

for s = 1:nScen
    name = scenarios{s}{1};
    idx_vec = scenarios{s}{2};          % state indices
    Ts_sec_vec = scenarios{s}{3};       % sampling times (sec), same length

    if numel(idx_vec) ~= numel(Ts_sec_vec)
        error('Scenario %d: idx_vec and Ts_sec_vec must have same length.', s);
    end

    % Compute multi-rate observability rank
    r = obs_rank_multirate(A, idx_vec, Ts_sec_vec, dt_sec, N_h);
    obs_ranks(s) = r;

    fprintf('\nScenario %2d: %s\n', s, name);
    fprintf('  Measured states: %s\n', mat2str(idx_vec));
    fprintf('  Sampling times (sec): %s\n', mat2str(Ts_sec_vec));
    fprintf('  Observability rank over %.2f h = %d (of %d)\n', ...
            T_horizon_hr, r, nx);
    if r == nx
        fprintf('  → FULLY OBSERVABLE (in this multi-rate sense).\n');
    else
        fprintf('  → NOT fully observable (some modes remain unobservable).\n');
    end
end

%% 8. Visualization: bar plot of observability rank per scenario
figure('Color','w','Name','Multi-rate observability ranks');
bar(obs_ranks, 'FaceColor',[0.2 0.4 0.8]);
hold on;
yline(nx,'r--','LineWidth',1.5);
hold off;

xlabel('Scenario index');
ylabel('Observability rank over horizon');
title(sprintf('Multi-rate observability over %.2f h (dt = %.1f s)', ...
    T_horizon_hr, dt_sec));
ylim([0 nx+0.5]);
grid on;

set(gca,'XTick',1:nScen);
names = cell(nScen,1);
for s = 1:nScen
    names{s} = scenarios{s}{1};
end
set(gca,'XTickLabel',names, 'XTickLabelRotation',30);

% Optional numeric labels on bars
for s = 1:nScen
    text(s, obs_ranks(s)+0.1, sprintf('%d', obs_ranks(s)), ...
        'HorizontalAlignment','center', 'Rotation',0, 'FontSize',8);
end

fprintf('\nVisualization generated: bar plot of multi-rate observability ranks.\n');

%% ===== Local function: multi-rate observability rank =====
function r = obs_rank_multirate(A, idx_vec, Ts_sec_vec, dt_sec, N_h)
% obs_rank_multirate
%   Computes observability rank for a continuous-time system dx/dt = A x
%   under multi-rate, discrete measurements of selected states.
%
%   idx_vec:       indices of measured states (e.g. [6 2 1])
%   Ts_sec_vec:    sampling periods of each measurement [sec]
%   dt_sec:        base time step [sec] for grid
%   N_h:           number of base steps over the horizon
%
%   The observability matrix is built as:
%     O = [ C_k0;
%           C_k1 * A_d^k1;
%           C_k2 * A_d^k2;
%           ... ]
%   where C_kj is the row selecting the measured state at sampling time k_j.

n  = size(A,1);
Ad = expm(A * (dt_sec/3600));   % discrete-time A for base step (dt in hours)

% Convert sampling periods to integer multiples of dt
Ns_vec = round(Ts_sec_vec / dt_sec);   % number of base steps between samples
m      = numel(idx_vec);

O = [];
Ak = eye(n);   % Ad^k

for k = 0:N_h
    % For each sensor, check if it samples at this time step
    for j = 1:m
        if mod(k, Ns_vec(j)) == 0
            % Build measurement row C_j: 1 in idx_vec(j), 0 elsewhere
            Cj = zeros(1,n);
            Cj(idx_vec(j)) = 1;
            % Append to observability matrix
            O = [O; Cj * Ak];
        end
    end
    % Advance Ak = Ad^(k+1)
    Ak = Ak * Ad;
end

r = rank(O);
end
