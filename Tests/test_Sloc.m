function  [SA, Sloc] = test_Sloc(A, params, mesh, matrices, B_mu, phi, p)

id     = ~mesh.id_dirichlet;
SA = zeros(mesh.npoint,1);

B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];

mu_fe = Valve_GetMu(B_ele,B_mu,params);
mu_inv = 1./((1-phi)*params.mu0 + (phi.^p).*mu_fe);
mu_inv = repmat(mu_inv,1,9);

Sloc  = sparse(matrices.ii(:),matrices.jj(:),(matrices.sloc_aa(:)).*mu_inv(:));

SA(id) = Sloc(id,id)*A(id);

end

