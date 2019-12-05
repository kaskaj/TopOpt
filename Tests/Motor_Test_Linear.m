%% Load data
clear all;
refin_level = 1;

folder_name = 'Motor_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrices');

model = [];
model.p         = 1;
model.nonlinear = 0;

rpm = 1500;
Tp = (pi/6)/((2*pi*rpm)/60);

time = (0:3*Tp/36:2*Tp);

for i = 1:length(time)
    
[J,phi] = Motor_MoveCurrent(mesh, params, i-1, time(i));

% % Plot prescribed domains
% figure;
% Motor_PlotEdges(params, 1, 'k-', 2);
% plot(mesh.x_mid(phi==0),mesh.y_mid(phi==0),'o','MarkerFaceColor','b');
% plot(mesh.x_mid(phi==1),mesh.y_mid(phi==1),'ro','MarkerFaceColor','r');
% axis equal;
 
% PlotData(mesh.x,mesh.y,mesh.elems2nodes,J);
% Motor_PlotEdges(params,max(J));

[A, B] = Motor_GetJ(phi, J,  mesh, matrices, params, model);

PlotData(mesh.x,mesh.y,mesh.elems2nodes,A);
caxis([-5e-4,5e-4]);
Motor_PlotEdges(params,max(A));

end


