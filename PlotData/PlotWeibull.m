function [] = PlotWeibull(model,params)

a = model.B_mu.a_w;

x = linspace(0,3,1000);
y = params.mu0 + a(3).*(x-a(1)).^(a(2)-1).*exp(-(x-a(1)).^a(2));
plot(x,y);
hold on;
axis([0,x(end),0,2e-3]);
xlabel('B^2');
ylabel('\mu');
grid on;
end

