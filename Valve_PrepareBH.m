% clear;
% load('Valve_Data/BHcurve.mat');
% B = BHcurve(:,2).^2;
% H = BHcurve(:,1);
% 
% mu = B./H;
% mu(1) = max(mu);        %Avoid NaN
% 
% f = fit(B,mu,'cubicinterp');
% 
% B_mu.max_mu = max(mu);
% B_mu.min_mu = min(mu);
% B_mu.max_B = B(2);
% B_mu.min_B = B(end);
% B_mu.f = f;
% 
% file_name   = fullfile('Valve_Data', 'B_mu.mat');
% save(file_name,'B_mu');


clear;
load('Valve_Data/BHcurve.mat');
B = BHcurve(:,2);
H = BHcurve(:,1);

mu = B./H;
mu(1) = max(mu);        %Avoid NaN

f = fit(B.^2,mu,'cubicinterp');

B_mu.max_mu = max(mu);
B_mu.min_mu = min(mu);
B_mu.max_B = B(2)^2;
B_mu.min_B = B(end)^2;
B_mu.f = f;

file_name   = fullfile('Valve_Data', 'B_mu.mat');
save(file_name,'B_mu');
