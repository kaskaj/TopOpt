%% Load data

refin_level = 4;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, 'B_mu_weib'),'B_mu_weib');
load(fullfile(folder_name, 'B_mu_exp'),'B_mu_exp');
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

%% Find the force

model = [];
model.p         = 1;
model.coil      = 1;
model.nonlinear = 1;
model.B_mu      = B_mu_weib;

[F, A, B, ~, ~, ~] = Valve_GetJ(phi, mesh, matrices, params, model);

fprintf('Force for nonlinear model: Fy = %d\n',F);

% PlotData(mesh.x,mesh.y,mesh.elems2nodes,A);
% Valve_PlotEdges(params,max(A));


