function [SA, Sloc] = Valve_GetSA(A, params, matrices, B_mu, phi, p)

B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];
mu_fe = Valve_GetMu(B_ele,B_mu,params);  

mu_inv = 1./((1-phi)*params.mu0 + (phi.^p).*mu_fe);
mu_inv = repmat(mu_inv,1,9);

Sloc  = sparse(matrices.ii(:),matrices.jj(:),(matrices.sloc_aa(:)).*mu_inv(:));

SA = Sloc*A;

end

