function y = h_meas(x,u,p)
% Available keys (order independent):
% 'h','C_DO','C_X_proxy','V','C_X','C_LAC','C_N','C_MEV'
if ~isfield(p,'meas_outputs'), error('Set p.meas_outputs as cellstr'); end
keys = p.meas_outputs; y = zeros(numel(keys),1);

V = x(1); C_X = x(2); C_LAC = x(3); C_N = x(4); C_MEV = x(5); C_DO = x(6);
h = (V/1000) / p.A_xsec;          % m
C_X_proxy = p.gamma0 + p.gamma1*C_X;

for i=1:numel(keys)
    switch keys{i}
        case 'h',          y(i)=h;
        case 'C_DO',       y(i)=C_DO;
        case 'C_X_proxy',  y(i)=C_X_proxy;
        case 'V',          y(i)=V;
        case 'C_X',        y(i)=C_X;
        case 'C_LAC',      y(i)=C_LAC;
        case 'C_N',        y(i)=C_N;
        case 'C_MEV',      y(i)=C_MEV;
        otherwise, error('Unknown meas key: %s', keys{i});
    end
end
end
