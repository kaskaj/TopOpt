%% Load data

refin_level = 1;

folder_name = 'Motor_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrices');

%% Prescribe fixed Air and Iron domains

phi = zeros(mesh.nelement,1);

ii_fix0 = ismember(mesh.tnum, [4,6]);
ii_fix1 = ismember(mesh.tnum, [1,2,3,5]);

%Coils
for i = 7:6:72    
    ii_fix0 = ii_fix0 | ismember(mesh.tnum, i:i+2);
end

%Iron
for i = 10:6:72    
    ii_fix1 = ii_fix1 | ismember(mesh.tnum, i:i+2);
end

ii_fix  = ii_fix0 | ii_fix1;
ii_opt  = ~ii_fix;

phi(ii_fix0)  = 0;
phi(ii_fix1)  = 1;

%Plot prescribed domains
Motor_PlotEdges(params, 1, 'k-', 2);
plot(mesh.x_mid(ii_fix0),mesh.y_mid(ii_fix0),'o','MarkerFaceColor','b');
plot(mesh.x_mid(ii_fix1),mesh.y_mid(ii_fix1),'ro','MarkerFaceColor','r');
axis equal;

figure;
plot(mesh.x(find(matrices.J)),mesh.y(find(matrices.J)),'o');

%Plot boundaries
% figure;
% plot(mesh.x(mesh.id_dirichlet==1),mesh.y(mesh.id_dirichlet==1),'o');
% hold on;
% plot(mesh.x(mesh.id_s1==1),mesh.y(mesh.id_s1==1),'ro');
% plot(mesh.x(mesh.id_s2==1),mesh.y(mesh.id_s2==1),'go');
% axis equal;

model = [];
model.p         = 1;
model.nonlinear = 0;

[A, B] = Motor_GetJ(phi, mesh, matrices, params, model);

PlotData(mesh.x,mesh.y,mesh.elems2nodes,A);
Motor_PlotEdges(params,max(A));


% normB = sqrt(B(:,1).^2 + B(:,2).^2);
% 
% PlotData(mesh.x,mesh.y,mesh.elems2nodes,normB);
% Motor_PlotEdges(params,max(normB));

