function dP = Valve_GetdP(phi, A, B, B_ele, Sloc_mu, mu_fe, dmu_fe, mesh, matrices, params, model)

id     = ~mesh.id_dirichlet;
npoint = mesh.npoint;
B_mu = model.B_mu;

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

dSA  = Valve_GetdSA(phi, A, B_ele, mesh, matrices, params, model, mu_fe, dmu_fe);
dSAf = Sloc_mu(id,id) + dSA(id,id);

alpha(id) = dSAf\f(id);


%% Compute derivative

B0    = (B_ele(:,1).^2) + (B_ele(:,2).^2);

da1 = (phi.^p).*(B_mu.a_w(3).*exp(-(B0-B_mu.a_w(1)).^B_mu.a_w(2)).*(B_mu.a_w(2).*(B0-B_mu.a_w(1)).^(2*B_mu.a_w(2)-2) - (B_mu.a_w(2)-1).*(B0-B_mu.a_w(1)).^(B_mu.a_w(2)-2)))./((1-phi)*mu_air + (phi.^p).*mu_fe).^2;
da2 = (phi.^p).*(B_mu.a_w(3).*exp(-(B0-B_mu.a_w(1)).^B_mu.a_w(2)).*(log(B0-B_mu.a_w(1)).*(B0-B_mu.a_w(1)).^(B_mu.a_w(2)-1) - log(B0-B_mu.a_w(1)).*(B0-B_mu.a_w(1)).^(2*B_mu.a_w(2)-1)))./((1-phi)*mu_air + (phi.^p).*mu_fe).^2;
da3 = (phi.^p).*((B0-B_mu.a_w(1)).^(B_mu.a_w(2)-1).*exp(-(B0-B_mu.a_w(1)).^B_mu.a_w(2)))./((1-phi)*mu_air + (phi.^p).*mu_fe).^2;

da1 = reshape(repmat(da1,1,9)', [], 1);
da2 = reshape(repmat(da2,1,9)', [], 1);
da3 = reshape(repmat(da3,1,9)', [], 1);

x1 = A(mesh.elems2nodes);
x2 = alpha(mesh.elems2nodes);
y1 = reshape(repmat(x1, 1, 3)', [], 1);
y2 = reshape([repmat(x2(:,1), 1, 3), repmat(x2(:,2), 1, 3), repmat(x2(:,3), 1, 3)]', [], 1);

dP1 = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* da1;
dP2 = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* da2;
dP3 = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* da3;

dP1 = sum(reshape(dP1, 9, []))';
dP2 = sum(reshape(dP2, 9, []))';
dP3 = sum(reshape(dP3, 9, []))';

dP = [dP1,dP2,dP3];

end

