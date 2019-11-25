load(fullfile(folder_name, 'B_mu'),'B_mu');
mu0 = 4*pi*1e-7;

a1 = B_mu.a_w(1);
a2 = B_mu.a_w(2);
a3 = B_mu.a_w(3);

a1_ = a1*0.8;
a2_ = a2*0.8;
a3_ = a3*0.8;


x = linspace(0,3,1000);
y0 = mu0 + a3.*(x-a1).^(a2-1).*exp(-(x-a1).^a2);
y1 = mu0 + a3.*(x-a1_).^(a2-1).*exp(-(x-a1_).^a2);
y2 = mu0 + a3.*(x-a1).^(a2_-1).*exp(-(x-a1).^a2_);
y3 = mu0 + a3_.*(x-a1).^(a2-1).*exp(-(x-a1).^a2);

plot(x,y0,'k-');
hold on;
plot(x,y1,'r-');
plot(x,y2,'g-');
plot(x,y3,'b-');
legend('Original','eps*a1','eps*a2','eps*a3');
grid on;