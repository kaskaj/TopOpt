clear all;

%% Load data
refin_level = 4;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, 'B_mu'),'B_mu');
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

p = 1;
coil = 1;
% n = 3;

mu_fe = params.mu0*params.mur*ones(mesh.nelement,1);

for i = 1:10
    
    [A, B_ele] = Valve_GetJ_nonlinear(phi, mesh, matrices, params, p, coil, mu_fe);
    mu_fe = Valve_GetMu(B_ele,B_mu);
    
    norm(A)
end

eroijoierjo

PlotData(mesh.x,mesh.y,mesh.elems2nodes,A);
Valve_PlotEdges(params,max(A));


% ele = delaunay(mesh.x_mid,mesh.y_mid);
% 
% PlotData(mesh.x,mesh.y,mesh.elems2nodes,B(:,1))
% caxis([min(B_ele(:,1)), max(B_ele(:,1))]);
% Valve_PlotEdges(params,max(B(:,1)))
% 
% PlotData(mesh.x_mid,mesh.y_mid,ele,B_ele(:,1))
% caxis([min(B_ele(:,1)), max(B_ele(:,1))]);
% Valve_PlotEdges(params,max(B_ele(:,1)))

