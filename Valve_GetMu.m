function [mu] = Valve_GetMu(B,B_mu)

B = (B(:,1).^2) + (B(:,2).^2);

ii_linear1 = (B < B_mu.max_B);
ii_linear2 = (B > B_mu.min_B);

ii_poly = ~ii_linear1 & ~ii_linear2;

mu = zeros(length(B),1);
mu(ii_linear1) = B_mu.max_mu;
mu(ii_linear2) = B_mu.min_mu;

mu(ii_poly) = B_mu.f(B(ii_poly));



end

