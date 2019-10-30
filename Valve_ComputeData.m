clear all;
close all;

addpath(genpath('./3rd_party'));

%% Load data

refin_level = 7;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)));

%% Extract parameters

t   = params.t;
mu0 = params.mu0;

x           = mesh.x;
y           = mesh.y;
id          = ~mesh.id_dirichlet;
npoint      = mesh.npoint;
elems2nodes = mesh.elems2nodes;

%% Prescribe current density

J = zeros(npoint,1);

J_nodes_1 = (x >= (params.x_c_min + 2*t)) & (x <= (params.x_c_max - t)) & (y >= (params.y_c_min + t)) & (y <= (params.y_c_max - t));
J_nodes_2 = (x >= (params.x_c_min + 2*t)) & (x <= (params.x_c_max - t)) & (y >= -(params.y_c_max - t)) & (y <= -(params.y_c_min + t));

J(J_nodes_1) = params.J_coil1;
J(J_nodes_2) = params.J_coil2;

%% Solve the system

f     = Mloc*J;
A     = zeros(npoint,1);
A(id) = Sloc_mu(id,id) \ f(id);
B     = [Mloc\(Clocy*A),-Mloc\(Clocx*A)];

%% Plot fields

field = A;
PlotData(x, y, elems2nodes, field);
Valve_PlotEdges(params, max(field));

field = B(:,1);
PlotData(x, y, elems2nodes, field);
Valve_PlotEdges(params, max(field));

%% Compute the force acting on the plunger (way 1)

F_x_aux = B(:,1)'*Clocx_plunger*B(:,1) - B(:,2)'*Clocx_plunger*B(:,2) + B(:,2)'*Clocy_plunger*B(:,1) + B(:,1)'*Clocy_plunger*B(:,2);
F_y_aux = -B(:,1)'*Clocy_plunger*B(:,1) + B(:,2)'*Clocy_plunger*B(:,2) + B(:,2)'*Clocx_plunger*B(:,1) + B(:,1)'*Clocx_plunger*B(:,2);
F       = 1/mu0*[F_x_aux; F_y_aux];

%% Compute the force acting on the plunger (way 2)

F_mod = zeros(4,2);
for i=1:4
    if i==1
        ii = y == params.y_piston_min & x >= params.x_piston_min & x <= params.x_piston_max;
    elseif i==2
        ii = x == params.x_piston_max & y >= params.y_piston_min & y <= params.y_piston_max;
    elseif i==3
        ii = y == params.y_piston_max & x >= params.x_piston_min & x <= params.x_piston_max;
    else
        ii = x == params.x_piston_min & y >= params.y_piston_min & y <= params.y_piston_max;
    end
    x_p  = x(ii);
    y_p  = y(ii);
    Bx_p = B(ii,1);
    By_p = B(ii,2);
    
    if i==1 || i==3
        [~,jj] = sort(x_p);
    else
        [~,jj] = sort(y_p);
    end
    x_p    = x_p(jj);
    y_p    = y_p(jj);
    Bx_p   = Bx_p(jj);
    By_p   = By_p(jj);
    
    T_p  = {1/2*(Bx_p.*Bx_p-By_p.*By_p), Bx_p.*By_p; Bx_p.*By_p, 1/2*(By_p.*By_p-Bx_p.*Bx_p)};
    if i==1
        T_pn = [-T_p{1,2}, -T_p{2,2}];
    elseif i==2
        T_pn = [T_p{1,1}, T_p{2,1}];
    elseif i==3
        T_pn = [T_p{1,2}, T_p{2,2}];
    else
        T_pn = [-T_p{1,1}, -T_p{2,1}];
    end
    xy_diff = diff([x_p, y_p]);
    xy_norm = sqrt(xy_diff(:,1).^2 + xy_diff(:,2).^2);
    
    T_mean  = 0.5*(T_pn(1:end-1,:) + T_pn(2:end,:));
    
    F_mod(i,:) = sum(xy_norm.*T_mean)';
end
F_mod = F_mod / params.mu0;
F_mod = sum(F_mod)';


F
F_mod



