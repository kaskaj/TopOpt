n = 100000;
H_new = linspace(0,max(H),n);
B_new = interp1(H,B,H_new);

mu = B_new./H_new;

i = 210;
p = polyfit(B_new(i:end),mu(i:end),3);
mu2 = zeros(n,1);

mu2(1:i) = max(mu);
mu2(i:end) = p(1)*B_new(i:end).^3 + p(2)*B_new(i:end).^2 + p(3)*B_new(i:end) + p(4);

figure,
%plot(B_new,mu);
plot(B_new,mu,B_new,mu2);
grid on
max_mu = max(mu);


B_mu.max_mu = max_mu;
B_mu.p = p;


file_name   = fullfile('Valve_Data', 'B_mu.mat');
save(file_name,'B_mu');
