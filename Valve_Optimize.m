clear all;

%% Load data

refin_level = 4;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, 'B_mu_exp'),'B_mu_exp');
load(fullfile(folder_name, 'B_mu_weib'),'B_mu_weib');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrices');

steps  = 30;
lambda = 1e-2;

model = [];
model.p         = 1;
model.coil      = 1;
model.nonlinear = 0;
model.B_mu      = B_mu_weib;

%% Prescribe fixed Air and Iron domains

phi = zeros(mesh.nelement,steps);

ii_fix0 = ismember(mesh.tnum, [2,3,5,6,7,8,9,10,14,15]);
ii_fix1 = ismember(mesh.tnum, 4);
ii_fix  = ii_fix0 | ii_fix1;
ii_opt  = ~ii_fix;

%Initial Iron
ii_init = (mesh.x_mid >= params.x_init_min) & (mesh.x_mid <= params.x_init_max) & (mesh.y_mid >= params.y_init_min) & (mesh.y_mid <= params.y_init_max);
phi(ii_init, 1) = 1;
phi(ii_fix0, 1)  = 0;
phi(ii_fix1, 1)  = 1;

%Plot prescribed domains
% figure;
% Valve_PlotEdges(params, 1);
% plot(mesh.x_mid(ii_fix0,1),mesh.y_mid(ii_fix0,1),'o','MarkerFaceColor','b');
% plot(mesh.x_mid(ii_fix1,1),mesh.y_mid(ii_fix1,1),'ro','MarkerFaceColor','r');
% plot(mesh.x_mid(ii_init,1),mesh.y_mid(ii_init,1),'ro','MarkerFaceColor','r');
% axis equal;

%% Optimization LOOP

F = zeros(steps,1);
F_round = nan(steps,1);
A = [];

dir = zeros(1,3);
eps = [1e-2,1e-2,1e-4];
model_new = model;
der = zeros(steps,3);
for i = 1:steps
    
    [F(i), A, B, B_ele, Sloc_mu, mu_fe, dmu_fe, f] = Valve_GetJ(phi(:,i), mesh, matrices, params, model_new, A);
    [dJdphi,dJdp] = Valve_GetdJ(phi(:,i), A, B, B_ele, Sloc_mu, mu_fe, dmu_fe, mesh, matrices, params, model_new);
    
    der(i,:) = dJdp;
    
    %Modify B-mu parameters
    a = model.B_mu.a_w;    
    dir(dJdp > 0) = 1;
    dir(dJdp < 0) = -1;
    a = a + dir.*eps;
    model_new.B_mu.a_w = a;

    
    phi(ii_fix0,i+1) = 0;
    phi(ii_fix1,i+1) = 1;
    phi(ii_opt,i+1)  = phi(ii_opt,i) - lambda*dJdphi(ii_opt);
    phi(ii_opt,i+1)  = max(min(phi(ii_opt,i+1), 1), 0);
    
    if mod(i, 10) == 0
        F_round(i) = Valve_GetJ(round(phi(:,i+1)), mesh, matrices, params, model_new);
        fprintf('step%d: Fy = %d, Fy_round = %d.\n',i, F(i), F_round(i));
    else
        fprintf('step%d: Fy = %d\n',i,F(i));
    end
end

model_true   = model;
model_true.p = 1;
model.B_mu.a_w = model.B_mu.a_w;

phi_final = round(phi(:,end));

F1 = Valve_GetJ(phi(:,end), mesh, matrices, params, model_true);
[F2, A, B] = Valve_GetJ(phi_final, mesh, matrices, params, model_true);

fprintf('Last step: Fy = %d, Fy_round = %d.\n', F1, F2);

ii1 = phi_final == 0;
ii2 = phi_final == 1;
ii3 = ~(ii1 | ii2);

%Plot results
SubPlotPhi(mesh, params, phi(:,end));

PlotData(mesh.x,mesh.y,mesh.elems2nodes,A)
Valve_PlotEdges(params,max(A))

normB = sqrt(B(:,1).^2 + B(:,2).^2);
PlotData(mesh.x,mesh.y,mesh.elems2nodes,normB);
Valve_PlotEdges(params,max(normB));
quiver3(mesh.x(1:50:end), mesh.y(1:50:end),max(normB)*ones(length(mesh.x(1:50:end)),1),B(1:50:end,1),B(1:50:end,2),zeros(length(mesh.x(1:50:end)),1),'white');

% Save results
% file_name = fullfile('Valve_Results', 'Fy2.mat');
% save(file_name, 'F', 'F_round', 'phi');




