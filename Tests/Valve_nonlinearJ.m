function [A, B, B_ele, mu_fe, dmu_fe, Sloc_mu] = Valve_nonlinearJ(phi, mesh, matrices, params, B_mu, p, coil, A0)

nonlinear = 1;
maxsteps = 100;
%Initial guess for Newton
mu_fe = params.mu0*params.mur*ones(mesh.nelement,1);

% if nargin < 7 || isempty(A0)
    [~, A, ~, B_ele, Sloc_mu, f] = Valve_GetJ(phi, mesh, matrices, params, p, coil, nonlinear, mu_fe);
% else

% Newton
[mu_fe, dmu_fe] = Valve_Newton(A, B_ele, phi, p, f, mesh, matrices, params, B_mu, maxsteps);
%Get A,B
[~, A, B, B_ele] = Valve_GetJ(phi, mesh, matrices, params, p, coil, nonlinear, mu_fe);

end

