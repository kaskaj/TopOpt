clear all;

load('Valve_Data/BHcurve.mat');
load(fullfile('Valve_Data', 'Param'), 'params');
B = BHcurve(:,2);
H = BHcurve(:,1);

mu0 = params.mu0;
mu = B./H;
mu(1) = max(mu);        %Avoid NaN

a = 3.8;

%Plot:
% mu2 = mu0 + (max(mu) - mu0)*exp(-a*(B-1));
% plot(B,mu);
% hold on;
% plot(B(2:end),mu2(2:end),'r-');
% grid on;

B_mu.max_mu = max(mu);
B_mu.min_mu = min(mu);
B_mu.max_B = B(2)^2;
B_mu.min_B = B(end)^2;
B_mu.a = a;

file_name   = fullfile('Valve_Data', 'B_mu.mat');
save(file_name,'B_mu');
