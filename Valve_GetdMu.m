function dMu = Valve_GetdMu(B_ele,B_mu,params)

B = (B_ele(:,1).^2) + (B_ele(:,2).^2);

ii_linear = (B < B_mu.max_B);
ii_exp = ~ii_linear;

dMu = zeros(length(B),1);

dMu(ii_linear) = 0;
dMu(ii_exp) = -B_mu.a*(B_mu.max_mu - params.mu0)*exp(-B_mu.a*(B(ii_exp)-B_mu.max_B));

end

