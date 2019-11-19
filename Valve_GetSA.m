function [SA, Sloc, mu_fe] = Valve_GetSA(phi, A, params, matrices, model)

p    = model.p;
B_mu = model.B_mu;

B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];
mu_fe = Valve_GetMu(B_ele,B_mu,params);  

mu_inv = 1./((1-phi)*params.mu0 + (phi.^p).*mu_fe);
mu_inv = repmat(mu_inv,1,9);

Sloc  = sparse(matrices.ii(:),matrices.jj(:),(matrices.sloc_aa(:)).*mu_inv(:));

SA = Sloc*A;

end

