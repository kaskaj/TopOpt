clear all;

%% Load data

refin_level = 4;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrices');

%% Prescribe fixed Air and Iron domains

phi = zeros(mesh.nelement,1);

ii_fix0 = ismember(mesh.tnum, [2,3,5,6,7,8,9,10,14,15]);
ii_fix1 = ismember(mesh.tnum, [4,11,12,13]);
ii_fix  = ii_fix0 | ii_fix1;
ii_opt  = ~ii_fix;

%Initial Iron
%ii_init = (mesh.x_mid >= params.x_init_min) & (mesh.x_mid <= params.x_init_max) & (mesh.y_mid >= params.y_init_min) & (mesh.y_mid <= params.y_init_max);
%phi(ii_init) = 1;
phi(ii_fix0)  = 0;
phi(ii_fix1)  = 1;

%Plot prescribed domains
figure;
Valve_PlotEdges(params, 1);
plot(mesh.x_mid(ii_fix0),mesh.y_mid(ii_fix0),'o','MarkerFaceColor','b');
plot(mesh.x_mid(ii_fix1),mesh.y_mid(ii_fix1),'ro','MarkerFaceColor','r');
%plot(mesh.x_mid(ii_init),mesh.y_mid(ii_init),'ro','MarkerFaceColor','r');
axis equal;

figure;
Valve_PlotEdges(params, 1);
plot(mesh.x(matrices.Br > 0),mesh.y(matrices.Br > 0),'o','MarkerFaceColor','b');
axis equal;

p = 1;
coil = 0;
[F, A, B, Sloc_mu] = Valve_GetJ(phi, mesh, matrices, params, p , coil);

PlotData(mesh.x,mesh.y,mesh.elems2nodes,A)
Valve_PlotEdges(params,max(A))

