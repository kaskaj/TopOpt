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

phi(ii_fix0)  = 0;
phi(ii_fix1)  = 1;

%Plot prescribed domains
% figure;
% Valve_PlotEdges(params, 1);
% plot(mesh.x_mid(ii_fix0),mesh.y_mid(ii_fix0),'o','MarkerFaceColor','b');
% plot(mesh.x_mid(ii_fix1),mesh.y_mid(ii_fix1),'ro','MarkerFaceColor','r');
% axis equal;
% 
% figure;
% Valve_PlotEdges(params, 1);
% plot(mesh.x(matrices.Br > 0),mesh.y(matrices.Br > 0),'o','MarkerFaceColor','b');
% axis equal;

model = [];
model.p         = 1;
model.coil      = 1;
model.nonlinear = 0;

[F, A, B, B_ele, Sloc_mu] = Valve_GetJ(phi, mesh, matrices, params, model);

fprintf('Force for linear model: Fy = %d\n',F);

PlotData(mesh.x,mesh.y,mesh.elems2nodes,A);
Valve_PlotEdges(params,max(A));

% ele = delaunay(mesh.x_mid,mesh.y_mid);
% 
% PlotData(mesh.x_mid,mesh.y_mid,ele,B_ele(:,1))
% caxis([min(B_ele(:,1)), max(B_ele(:,1))]);
% Valve_PlotEdges(params,max(B_ele(:,1)))

