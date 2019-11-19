function [mu, dMu] = Valve_GetMu(B_ele, params, model)

B_mu = model.B_mu;
B    = (B_ele(:,1).^2) + (B_ele(:,2).^2);

ii_linear = (B < B_mu.max_B);
ii_exp = ~ii_linear;

mu  = zeros(length(B),1);
dMu = zeros(length(B),1);

mu(ii_linear)  = B_mu.max_mu;
mu(ii_exp)     = params.mu0 + (B_mu.max_mu - params.mu0)*exp(-B_mu.a*(B(ii_exp)-B_mu.max_B));
dMu(ii_linear) = 0;
dMu(ii_exp)    = -B_mu.a*(B_mu.max_mu - params.mu0)*exp(-B_mu.a*(B(ii_exp)-B_mu.max_B));

end

