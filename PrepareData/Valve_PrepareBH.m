clear all;

load('Valve_Data/BHcurve.mat');
load(fullfile('Valve_Data', 'Param'), 'params');
B = BHcurve(:,2);
H = BHcurve(:,1);

mu0 = params.mu0;
mu = B./H;
mu(1) = max(mu);        %Avoid NaN

%% Poly:
expo = 3.8;

%% Weibull
Ax = 0;     Ay = 1.4e-3;
Bx = 0.5;	By = 1.5e-3;
Cx = 1;     Cy = 1.3e-3;

x = [Ax, Bx, Cx, B(3:end)'];
y = [Ay, By, Cy, mu(3:end)'];

a1 = -10.^linspace(0,2,40);
a2 = 10.^linspace(1,3,40);
a3 = 10.^linspace(1,6,40);
a4 = 10.^linspace(-4,0,40);

a = combvec(a1,a2,a3,a4);
a = a';

pred = zeros(length(x),size(a,1));
for i=1:length(x)
pred(i,:) = mu0 + a(:,4).*(((x(i)-a(:,1))./a(:,2)).^(a(:,3)-1)).*exp(-(((x(i)-a(:,1))./a(:,2)).^a(:,3)));
end

solv = pred' - repmat(y, size(a,1), 1);
[m,ii] = min(sum(abs([5*solv(:,1), 9*solv(:,2:3), solv(:,4:end)]),2));

fprintf('Fitting error (Weibull) = %d\n',m);

weib = a(ii,:);

%% Plot

%Plot:
B_long = linspace(0,B(end),1000);
mu2 = mu0 + (max(mu) - mu0)*exp(-expo*(B-1));
mu3 = mu0 + weib(4).*(((B_long-weib(1))./weib(2)).^(weib(3)-1)).*exp(-(((B_long-weib(1))./weib(2)).^weib(3)));
plot(B,mu,'k-');
hold on;
plot(B(2:end),mu2(2:end),'r-');
hold on;
plot(B_long,mu3,'b-');
axis([0,B(end),0,2e-3]);
legend('Comsol','Exponential','Weibull');
grid on;

B_mu.max_mu = max(mu);
B_mu.min_mu = min(mu);
B_mu.max_B = B(2)^2;
B_mu.min_B = B(end)^2;
B_mu.a_e = expo;
B_mu.a_w = weib;

file_name   = fullfile('Valve_Data', 'B_mu.mat');
save(file_name,'B_mu');
