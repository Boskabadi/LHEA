function [A,B,C,D,f0,y0] = linearize_CSTRO(x,u,p)
% Linearization of CSTRO ODE
% ẋ = f(x,u)
% y  = h(x)

%% Evaluate model
f0 = CSTRO_ode(x,u,p);

nx = numel(x);
nu = numel(u);

A = zeros(nx,nx);
B = zeros(nx,nu);

%% Finite-difference steps (unit-scaled)
dx = 1e-5 * max(1, abs(x));    % preserves physical units
du = 1e-5 * max(1, abs(u));

%% A = ∂f/∂x
for i = 1:nx
    xp = x;
    xp(i) = xp(i) + dx(i);
    A(:,i) = (CSTRO_ode(xp,u,p) - f0) / dx(i);
end

%% B = ∂f/∂u
for j = 1:nu
    up = u;
    up(j) = up(j) + du(j);
    B(:,j) = (CSTRO_ode(x,u,p) - CSTRO_ode(x,up,p)) / du(j);
end

%% Measurement model
% y = [ h; C_DO; biomass proxy ]

V    = x(1);    % [L]
C_X  = x(2);    % [g/L]
C_DO = x(6);    % [g/L]

h = (V/1000) / p.A_xsec;   % [m]

y0 = [
    h;
    C_DO;
    p.gamma0 + p.gamma1*C_X
];

%% C = ∂h/∂x
C = zeros(3,nx);

% dh/dV [m/L]
C(1,1) = 1/(1000*p.A_xsec);

% d(C_DO)/dx
C(2,6) = 1;

% d(proxy)/dC_X
C(3,2) = p.gamma1;

%% D = ∂h/∂u
D = zeros(3,nu);
end
