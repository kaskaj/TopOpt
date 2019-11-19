function dJ = Valve_GetdJ(phi, A, B, B_ele, Sloc_mu, mu_fe, dmu_fe, mesh, matrices, params, model)

id     = ~mesh.id_dirichlet;
npoint = mesh.npoint;

mu_air = params.mu0;
p      = model.p;

%% Compute beta and gamma

Cp_x  = (-matrices.Clocy_plunger' - matrices.Clocy_plunger)*B(:,1) + (matrices.Clocx_plunger' + matrices.Clocx_plunger)*B(:,2);
Cp_y  = (matrices.Clocx_plunger' + matrices.Clocx_plunger)*B(:,1) + (matrices.Clocy_plunger' + matrices.Clocy_plunger)*B(:,2);
beta  = -(1/params.mu0) * (matrices.Mloc\Cp_x);
gamma = -(1/params.mu0) * (matrices.Mloc\Cp_y);

%% Compute alpha

f         = matrices.Clocy'*beta - matrices.Clocx'*gamma;
alpha     = zeros(npoint,1);

if model.nonlinear == 1
    dSA  = Valve_GetdSA(phi, A, B_ele, mesh, matrices, params, model, mu_fe, dmu_fe);
    dSAf = Sloc_mu(id,id) + dSA(id,id);
    
    alpha(id) = dSAf\f(id);    
else
    alpha(id) = Sloc_mu(id,id)\f(id);
end

%% Compute derivative

dmu_inv = repmat((mu_air - p.*mu_fe.*phi.^(p-1))./((1-phi)*mu_air + (phi.^p).*mu_fe).^2, 1, 9);
dmu_inv = reshape(dmu_inv', [], 1);

x1 = A(mesh.elems2nodes);
x2 = alpha(mesh.elems2nodes);
y1 = reshape(repmat(x1, 1, 3)', [], 1);
y2 = reshape([repmat(x2(:,1), 1, 3), repmat(x2(:,2), 1, 3), repmat(x2(:,3), 1, 3)]', [], 1);

dJ = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* dmu_inv;
dJ = sum(reshape(dJ, 9, []))';
dJ = -dJ;

end

