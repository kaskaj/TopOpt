clear all;

%% Load data
refin_level = 6;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrices');

steps = 3;
p = 1.0;
lambda = 1e-1;

%% Prescribe fixed Air and Iron domains

phi = zeros(mesh.nelement,steps);

ii_fix0 = ismember(mesh.tnum, [2,3,5,6,7,8]);
ii_fix1 = ismember(mesh.tnum, 4);
ii_fix  = ii_fix0 | ii_fix1;
ii_opt  = ~ii_fix;

%Initial Iron
mu_nodes = (mesh.x_mid >= params.x_init_min) & (mesh.x_mid <= params.x_init_max) & (mesh.y_mid >= params.y_init_min) & (mesh.y_mid <= params.y_init_max);
phi(mu_nodes, 1) = 1;
phi(ii_fix0, 1)  = 0;
phi(ii_fix1, 1)  = 1;

%Plot initial Iron
% figure;
% Valve_PlotEdges(params, 1);
% plot(mesh.x_mid(ii_fix,1),mesh.y_mid(ii_fix,1),'ko','MarkerFaceColor','k');
% axis equal;

%% Optimization LOOP

F = zeros(steps,1);
F_round = nan(steps,1);

for i = 1:steps
    [F(i), A, B, Sloc_mu] = Valve_GetJ(phi(:,i), mesh, matrices, params, p);
    dJ = Valve_GetdJ(phi(:,i), Sloc_mu, A, B, mesh, matrices, params, p);
    
    phi(ii_fix0,i+1) = 0;
    phi(ii_fix1,i+1) = 1;
    phi(ii_opt,i+1)  = phi(ii_opt,i) - lambda*dJ(ii_opt);
    phi(ii_opt,i+1)  = max(min(phi(ii_opt,i+1), 1), 0);
    
    if mod(i, 10) == 0
        F_round(i) = Valve_GetJ(round(phi(:,i+1)), mesh, matrices, params, p);
        fprintf('step%d: Fy = %d, Fy_round = %d.\n',i, F(i), F_round(i));
    else
        fprintf('step%d: Fy = %d\n',i,F(i));
    end
end

phi_final = round(phi(:,end));

F1 = Valve_GetJ(phi(:,end), mesh, matrices, params, 1);
[F2, A, B] = Valve_GetJ(phi_final, mesh, matrices, params, 1);

fprintf('Last step: Fy = %d, Fy_round = %d.\n', F1, F2);

ii1 = phi_final == 0;
ii2 = phi_final == 1;
ii3 = ~(ii1 | ii2);

%Plot results
SubPlotPhi(mesh, params, phi(:,end));

%Save results
%file_name = fullfile('Valve_Results', 'Fy.mat');
%save(file_name, 'F', 'F_round', 'phi');




