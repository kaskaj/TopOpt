%% Load data

refin_level = 2;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, 'B_mu'),'B_mu');
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
model.B_mu      = B_mu;
model.aprox     = 'Weibull';

%% Compute derivatives

[F, A0, B, B_ele, Sloc_mu, mu_fe, dmu_fe] = Valve_GetJ(phi, mesh, matrices, params, model);
df_x = Valve_GetdJdp(phi, A0, B, B_ele, Sloc_mu, mu_fe, dmu_fe, mesh, matrices, params, model);


f = @(model) Valve_GetJ(phi, mesh, matrices, params, model);

for i = 1:3
    dfi = df_x(:,i);
    err = Diff_DerivativesStruct(f, dfi, model, i);

    fprintf('The relative error (dJdp)= %1.3e\n', err);
end



