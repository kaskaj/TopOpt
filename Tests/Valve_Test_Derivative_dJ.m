clear all;
close all;

%% Load data

refin_level = 4;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
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

p = 1;
coil = 1;
nonlinear = 0;

%% Compute derivatives

[F, A, B, ~, Sloc_mu] = Valve_GetJ(phi, mesh, matrices, params, p, coil, nonlinear);
dJ = Valve_GetdJ(phi, Sloc_mu, A, B, mesh, matrices, params, p, nonlinear);

f = @(phi) Valve_GetJ(phi, mesh, matrices, params, p, coil, nonlinear);
g = @(x,y) Valve_GetdJ(phi, Sloc_mu, A, B, mesh, matrices, params, p, nonlinear);

for i = 1:2
    if i == 1
        dir = dJ;
    else
        dir = sin(mesh.x_mid + mesh.y_mid);
    end
    err = Diff_Derivatives(f, g, phi, dir);

    fprintf('The relative error (dJ) = %1.3e\n', err);
end



