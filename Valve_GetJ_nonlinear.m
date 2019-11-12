function [A, B_ele, Sloc_mu] = Valve_GetJ_nonlinear(phi, mesh, matrices, params, p, coil, mu_fe)

id     = ~mesh.id_dirichlet;
npoint = mesh.npoint;

mu0 = params.mu0;

%% phi -> mu, dmu

mu_inv = 1./((1-phi)*mu0 + (phi.^p).*mu_fe);
mu_inv = repmat(mu_inv,1,9);

Sloc_mu  = sparse(matrices.ii(:),matrices.jj(:),(matrices.sloc_aa(:)).*mu_inv(:));

%% Solve the system
if coil == 1    %Turn on the current
    f = matrices.Mloc*matrices.J + (1/params.mu0)*matrices.Clocy*matrices.Br;
else            %Turn off the current
    f = (1/params.mu0)*matrices.Clocy*matrices.Br;
end

A     = zeros(npoint,1);
A(id) = Sloc_mu(id,id) \ f(id);
% B     = [matrices.Mloc\(matrices.Clocy*A),-matrices.Mloc\(matrices.Clocx*A)];
B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];

% F_x_aux = B(:,1)'*matrices.Clocx_plunger*B(:,1) - B(:,2)'*matrices.Clocx_plunger*B(:,2) + B(:,2)'*matrices.Clocy_plunger*B(:,1) + B(:,1)'*matrices.Clocy_plunger*B(:,2);
% F_y_aux = -B(:,1)'*matrices.Clocy_plunger*B(:,1) + B(:,2)'*matrices.Clocy_plunger*B(:,2) + B(:,2)'*matrices.Clocx_plunger*B(:,1) + B(:,1)'*matrices.Clocx_plunger*B(:,2);
% F       = 1/params.mu0*[F_x_aux; F_y_aux];

% F = F(2);   %force in y direction
% F = -F;

end

