function [mu] = Valve_GetMu(B,B_mu)

B = (B(:,1).^2) + (B(:,2).^2);

ii_linear = (B < 1);
ii_poly = ~ii_linear;

mu = zeros(length(B),1);
mu(ii_linear) = B_mu.max_mu;

mu(ii_poly) = B_mu.p(1)*B(ii_poly).^3 + B_mu.p(2)*B(ii_poly).^2 + B_mu.p(3)*B(ii_poly) + B_mu.p(4);


end

