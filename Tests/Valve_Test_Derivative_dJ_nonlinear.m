%% Load data

refin_level = 2;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, 'B_mu_exp'),'B_mu_exp');
load(fullfile(folder_name, 'B_mu_weib'),'B_mu_weib');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrices');

x_mid = mesh.x_mid;
y_mid = mesh.y_mid;

%% Prescribe fixed Air and Iron domains

phi = zeros(mesh.nelement,1);

ii_fix0 = ismember(mesh.tnum, [2,3,5,6,7,8,9,10,14,15]);
ii_fix1 = ismember(mesh.tnum, [4,11,12,13]);
ii_fix  = ii_fix0 | ii_fix1;
ii_opt  = ~ii_fix;

phi(ii_fix0)  = 0;
phi(ii_fix1)  = 1;

model = [];
model.p         = 1;
model.coil      = 1;
model.nonlinear = 1;
model.B_mu      = B_mu_weib;

%% Compute derivatives

[F, A, B, B_ele, Sloc_mu, mu_fe, dmu_fe] = Valve_GetJ(phi, mesh, matrices, params, model);
[df_x1,df_x2] = Valve_GetdJ(phi, A, B, B_ele, Sloc_mu, mu_fe, dmu_fe, mesh, matrices, params, model);

f = @(phi) Valve_GetJ(phi, mesh, matrices, params, model);

%dJdphi
for i = 1:2
    if i == 1
        dir = df_x1;
    else
        dir = sin(mesh.x_mid + mesh.y_mid);
    end
    err = Diff_Derivatives(f, df_x1, phi, dir);

    fprintf('The relative error (dJdphi) - nonlinear = %1.3e\n', err);
end

%dJdp
for i = 1:3
    f = @(a) PassStruct(a, i, phi, mesh, matrices, params, model);

    dfi = df_x2(:,i);
    err = Diff_Derivatives(f, dfi, model.B_mu.a_w(i), dfi);

    fprintf('The relative error (dJdp) %d = %1.3e\n',i , err);
end



