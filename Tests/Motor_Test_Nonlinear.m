%% Load data
clear all;
refin_level = 5;

folder_name = 'Motor_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, 'B_mu_weib'),'B_mu_weib');
load(fullfile(folder_name, 'B_mu_exp'),'B_mu_exp');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrices');

model = [];
model.p         = 1;
model.nonlinear = 1;
model.B_mu      = B_mu_weib;

rpm = 1500;
Tp = (pi/6)/((2*pi*rpm)/60);

time = (0:3*Tp/36:2*Tp);
Torque = zeros(length(time),4);

tic
for i = 1:length(time)
% for i = 1:1
fprintf('Step: %d of %d\n',i,length(time));      
[J,phi] = Motor_MoveCurrent(mesh, params, i-1, time(i));

% Plot prescribed domains
% figure;
% Motor_PlotEdges(params, 1, 'k-', 2);
% plot(mesh.x_mid(phi==0),mesh.y_mid(phi==0),'o','MarkerFaceColor','b');
% plot(mesh.x_mid(phi==1),mesh.y_mid(phi==1),'ro','MarkerFaceColor','r');
% axis equal;
 
% PlotData(mesh.x,mesh.y,mesh.elems2nodes,J);
% Motor_PlotEdges(params,max(J));

[A, B, T1, T2, T3, T4] = Motor_GetJ(phi, J,  mesh, matrices, params, model);
Torque(i,:) = [T1,T2,T3,T4];

% PlotData(mesh.x,mesh.y,mesh.elems2nodes,A);
% Motor_PlotEdges(params,max(A));

% normB = sqrt(B(:,1).^2 + B(:,2).^2);
% PlotData(mesh.x,mesh.y,mesh.elems2nodes,normB);
% Motor_PlotEdges(params,max(normB));

end
toc

figure
plot(time, Torque);
grid on;
