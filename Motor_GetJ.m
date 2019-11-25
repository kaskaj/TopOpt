function [A, B] = Motor_GetJ(phi, mesh, matrices, params, model)

id     = ~mesh.id_dirichlet;
npoint = mesh.npoint;

mu_air = params.mu0;
p      = model.p;

%% Right-hand side

f = matrices.Mloc*matrices.J;

%% Compute A

mu_fe   = params.mu0*params.mur;
mu_inv  = 1./((1-phi)*mu_air + (phi.^p).*mu_fe);
mu_inv  = repmat(mu_inv,1,9);
Sloc_mu = sparse(matrices.ii(:),matrices.jj(:),((matrices.sloc_aa(:)).*mu_inv(:)));

A     = zeros(npoint,1);
A(id) = Sloc_mu(id,id) \ f(id);


%% Compute B

B     = [matrices.Mloc\(matrices.Clocy*A),-matrices.Mloc\(matrices.Clocx*A)];

end

