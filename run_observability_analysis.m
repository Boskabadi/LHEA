%% run_observability_analysis.m
% Observability analysis for CSTRO under realistic measurement conditions
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
p.gamma0 = 0;                % offset for biomass proxy
p.gamma1 = 1;                % gain for biomass proxy (proxy ≈ C_X)

% Oxygen stoichiometry (example values; set from your stoichiometry if needed)
p.Y_O2_LAC_growth = 0.52;    % g O2 / g lactose for growth pathway
p.Y_O2_LAC_prod   = 0.12;    % g O2 / g lactose for product pathway

%% 3. Linearize the nonlinear model at (x0, u0)
% linearize_CSTRO returns:
%   A,B: state & input Jacobians
%   C_full,D: measurement Jacobians for [h; C_DO; biomass_proxy]
%   f0,y0: model and output at the linearization point
[A,B,C_full,D,f0,y0] = linearize_CSTRO(x0,u0,p);

nx = size(A,1);

%% 4. Define realistic measurement sets (“experiments”)
% State indices: [1 V, 2 C_X, 3 C_LAC, 4 C_N, 5 C_MEV, 6 C_DO]
%
% Each cell: { name, measured_state_indices }
measurement_sets = {
    {'Senario 1',                 [6]}                 % optical DO probe
    {'Senario 2',              [1]}                 % level sensor
    {'Senario 3',      [2]}                 % biomass via OD/capacitance
    {'Senario 4',              [6 1]}
    {'Senario 5',            [6 2]}
    {'Level + Biomass',         [1 2]}
    {'DO + Level + Biomass',    [6 1 2]}
    {'DO + Lactose',            [6 3]}              % DO + substrate soft-sensor
    {'Full soft-sensor set',    [1 2 3 6]}          % V, C_X, C_LAC, C_DO
    {'All measurable (no C_N)', [1 2 3 5 6]}        % V, C_X, C_LAC, C_MEV, C_DO
};

nExp = numel(measurement_sets);
obs_ranks = zeros(nExp,1);

%% 5. Loop through experiments and compute observability rank
fprintf('\n=== OBSERVABILITY ANALYSIS (linearized at x0,u0) ===\n');
for k = 1:nExp
    name = measurement_sets{k}{1};
    idx  = measurement_sets{k}{2};   % measured state indices

    % Build measurement matrix C_meas from direct state measurements
    C_meas = zeros(numel(idx), nx);
    for j = 1:numel(idx)
        C_meas(j, idx(j)) = 1;
    end

    % Observability rank
    r = obs_rank(A, C_meas);
    obs_ranks(k) = r;

    % Print result
    fprintf('\nExperiment %2d: %s\n', k, name);
    fprintf('  Measured states: %s\n', mat2str(idx));
    fprintf('  Observability rank = %d (of %d states)\n', r, nx);
    if r == nx
        fprintf('  → FULLY OBSERVABLE at this operating point.\n');
    else
        fprintf('  → NOT fully observable (some modes unobservable).\n');
    end
end

%% 6. Visualization: bar plot of observability rank per experiment
figure('Color','w','Name','Observability ranks across measurement sets');
bar(obs_ranks, 'FaceColor',[0.2 0.4 0.8]);
hold on;
yline(nx,'r--','LineWidth',1.5);
hold off;

xlabel('Experiment index');
ylabel('Observability rank');
title('Observability rank for different measurement configurations');
ylim([0 nx+0.5]);
grid on;

% Put experiment names on x-axis as tick labels (rotated)
set(gca,'XTick',1:nExp);
names = cell(nExp,1);
for k = 1:nExp
    names{k} = measurement_sets{k}{1};
end
set(gca,'XTickLabel',names, 'XTickLabelRotation',45);

%% 7. Optional: text labels on bars
for k = 1:nExp
    text(k, obs_ranks(k)+0.1, sprintf('%d', obs_ranks(k)), ...
        'HorizontalAlignment','center', 'Rotation',0, 'FontSize',8);
end

fprintf('\nVisualization generated: bar plot of observability ranks.\n');

%% ===== Local function for observability rank =====
function r = obs_rank(A, C)
% Compute observability matrix rank for (A,C)
% O = [ C; C*A; C*A^2; ...; C*A^(n-1) ]
n = size(A,1);
O = [];
Ak = eye(n);
for k = 0:n-1
    O = [O; C*Ak];
    Ak = Ak*A;
end
r = rank(O);
end
