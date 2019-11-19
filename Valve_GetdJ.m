function dJ = Valve_GetdJ(phi, Sloc_mu, A, B, B_ele, mesh, matrices, params, p, nonlinear, mu_fe, dmu_fe, B_mu)

id     = ~mesh.id_dirichlet;
npoint = mesh.npoint;

mu1 = params.mu0;

if nonlinear == 1
    mu2 = mu_fe;
else
    mu2 = params.mu0*params.mur;
end

%% Compute derivative

Cp_x  = (-matrices.Clocy_plunger' - matrices.Clocy_plunger)*B(:,1) + (matrices.Clocx_plunger' + matrices.Clocx_plunger)*B(:,2);
Cp_y  = (matrices.Clocx_plunger' + matrices.Clocx_plunger)*B(:,1) + (matrices.Clocy_plunger' + matrices.Clocy_plunger)*B(:,2);
beta  = -(1/params.mu0) * (matrices.Mloc\Cp_x);
gamma = -(1/params.mu0) * (matrices.Mloc\Cp_y);

f         = matrices.Clocy'*beta - matrices.Clocx'*gamma;
alpha     = zeros(npoint,1);

if nonlinear == 1
    dSA = Valve_GetdSA(A, B_ele, phi, mesh, matrices, params, B_mu, p, mu_fe, dmu_fe);
    dSAf = Sloc_mu(id,id) + dSA(id,id);
    
    alpha(id) = dSAf\f(id);    
else
    alpha(id) = Sloc_mu(id,id)\f(id);
end

dmu_inv = repmat((mu1 - p.*mu2.*phi.^(p-1))./((1-phi)*mu1 + (phi.^p).*mu2).^2, 1, 9);
dmu_inv = reshape(dmu_inv', [], 1);

x1 = A(mesh.elems2nodes);
x2 = alpha(mesh.elems2nodes);
y1 = reshape(repmat(x1, 1, 3)', [], 1);
y2 = reshape([repmat(x2(:,1), 1, 3), repmat(x2(:,2), 1, 3), repmat(x2(:,3), 1, 3)]', [], 1);

dJ = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* dmu_inv;
dJ = sum(reshape(dJ, 9, []))';
dJ = -dJ;

end

