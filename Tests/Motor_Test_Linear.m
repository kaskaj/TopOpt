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

time = linspace(0,Tp,37);

for i = 1:length(time)/2

[J,phi] = Motor_MoveCurrent(mesh, params, -i+1, time(i));
[A, B] = Motor_GetJ(phi, J,  mesh, matrices, params, model);

PlotData(mesh.x,mesh.y,mesh.elems2nodes,A);
caxis([-20e-4,20e-4]);
Motor_PlotEdges(params,max(A));

end


