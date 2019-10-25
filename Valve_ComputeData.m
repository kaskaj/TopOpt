clear all;

addpath(genpath('./3rd_party'));

%% Set parameters

levels = 3;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'param'));
% load(fullfile(folder_name, 'BHcurve'));
load(fullfile(folder_name, sprintf('mesh%d.mat', levels)));
load(fullfile(folder_name, sprintf('Matrices%d.mat', levels)));


% Prescribe current density

J = zeros(npoint,1);

J_nodes_1 = (x >= (x_c_min + 2*t)) & (x <= (x_c_max - t)) & (y >= (y_c_min + t)) & (y <= (y_c_max - t));
J_nodes_2 = (x >= (x_c_min + 2*t)) & (x <= (x_c_max - t)) & (y >= -(y_c_max - t)) & (y <= -(y_c_min + t));

J(J_nodes_1) =   J_coil1;
J(J_nodes_2) =   J_coil2;

A = zeros(npoint,1);
B = zeros(npoint,2);
normB = zeros(npoint,1);

f = Mloc*J;
A(id) = Sloc_mu(id,id) \ f(id);

B(:,:) = [-Mloc\(Clocy*A),Mloc\(Clocx*A)];
normB(:) = sqrt(B(:,1).^2 + B(:,2).^2);

PlotData(x,y,elems2nodes,A);
Valve_PlotEdges(max(A));
