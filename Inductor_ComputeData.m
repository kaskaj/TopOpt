clear all;

addpath(genpath('./3rd_party'));

%% Set parameters

dt = 1e-4;          %Time step
t_end = 2e-4;
levels = 5;

folder_name = 'Inductor_Data';

load(fullfile(folder_name, 'param'));
% load(fullfile(folder_name, 'BHcurve'));
load(fullfile(folder_name, sprintf('mesh%d.mat', levels)));
load(fullfile(folder_name, sprintf('Matrices%d.mat', levels)));


% Prescribe current density

J = zeros(npoint,1);

J_nodes_1 = (x >= (2*L - d2/2)) & (x <= (2*L + d2/2)) & (y >= (-d2/2)) & (y <= (d2/2));
J_nodes_2 = (x >= (-2*L - d2/2)) & (x <= (-2*L + d2/2)) & (y >= (-d2/2)) & (y <= (d2/2));

J(J_nodes_1) =  -J_coil;
J(J_nodes_2) =   J_coil;

%InductorInit();

time = linspace(0,t_end,t_end/dt);
npoints = length(x);

A_time = zeros(npoints,length(time));
J_time = zeros(npoints,length(time));

B_time = zeros(npoints,2,length(time));
normB_time = zeros(npoints,length(time));

step = 1;
for t = time
    
    %J_time(:,step) = J*cos(2*pi*f0*(t+dt/2));
    J_time(:,step) = J;
    
    A = zeros(npoints,1);
    %interp1(BHcurve(:,1),BHcurve(:,2),normB)
    
    if step == 1
        %f = Mloc*J_time(:,step);
        %A(id) = Sloc(id,id) \ f(id);
        %A_time(:,step) = A;
        A_time(:,step) = zeros(npoints,1);
    else
        S = (Sloc_mu + Mloc_gamma/dt);
        f = Mloc*J_time(:,step) + Mloc_gamma*A_time(:,step-1)/dt;
        A(id) = S(id,id) \ f(id);
        A_time(:,step) = A;
    end
    
    step = step + 1;
    
end

step = 2;

B_time(:,:,step) = [-Mloc\(Clocy*A_time(:,step)),Mloc\(Clocx*A_time(:,step))];
normB_time(:,step) = sqrt(B_time(:,1,step).^2 + B_time(:,2,step).^2);

PlotData(x,y,A_time(:,step));

PlotData(x,y,B_time(:,1,step));

%hold on;
%quiver(x,y,B_time(:,1,step),B_time(:,2,step),'color','white');


% x_bound = 600e-3;
% y_bound = 250e-3;
% y_line = 0;
% x_line = -75e-3;
% y_tol = 1e-4;
% x_tol = 1e-3;
%
% line = find(y < (y_line + y_tol) & y > (y_line - y_tol) & x > -x_bound & x < x_bound);
% line = find(x < (x_line + x_tol) & x > (x_line - x_tol) & y > -y_bound & y < y_bound);
% line = find(x < (x_line + x_tol) & x > (x_line - x_tol));
%
% A_line = A_time(line,step);
% figure;
% plot(nodes2coord(line,1),A_line,'o');
% grid on;
