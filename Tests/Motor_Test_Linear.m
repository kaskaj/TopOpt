%% Load data
clear all;
refin_level = 3;

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
Force_x1 = zeros(length(time),1);
Force_y1 = zeros(length(time),1);
Force_x2 = zeros(length(time),1);
Force_y2 = zeros(length(time),1);

for i = 1:length(time)
% for i = 1:1
    
[J,phi] = Motor_MoveCurrent(mesh, params, i-1, time(i));

% % Plot prescribed domains
% figure;
% Motor_PlotEdges(params, 1, 'k-', 2);
% plot(mesh.x_mid(phi==0),mesh.y_mid(phi==0),'o','MarkerFaceColor','b');
% plot(mesh.x_mid(phi==1),mesh.y_mid(phi==1),'ro','MarkerFaceColor','r');
% axis equal;
%  
% PlotData(mesh.x,mesh.y,mesh.elems2nodes,J);
% Motor_PlotEdges(params,max(J));

[A, B, Fx1, Fy1, Fx2, Fy2] = Motor_GetJ(phi, J,  mesh, matrices, params, model);
Force_x1(i) = Fx1;
Force_y1(i) = Fy1;
Force_x2(i) = Fx2;
Force_y2(i) = Fy2;

% PlotData(mesh.x,mesh.y,mesh.elems2nodes,A);
% caxis([-5e-4,5e-4]);
% Motor_PlotEdges(params,max(A));

% normB = sqrt(B(:,1).^2 + B(:,2).^2);
% PlotData(mesh.x,mesh.y,mesh.elems2nodes,normB);
% Motor_PlotEdges(params,max(normB));

end

figure
plot(time,Force_x1,time,Force_y1)

figure
plot(time,Force_x2,time,Force_y2)



