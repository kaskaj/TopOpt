clear all;
close all;

%% Load data

refin_level = 4;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, 'B_mu'), 'B_mu');
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
coil = 0;
nonlinear = 0;

%% Compute derivatives

[~, A0, ~, B_ele] = Valve_GetJ(phi, mesh, matrices, params, p, coil, nonlinear);

[SA, Sloc] = test_Sloc(A0, params,matrices,B_mu,phi,p);

f = @(A) test_Sloc(A, params,matrices,B_mu,phi,p) - Sloc*A;
g = @(A) Valve_GetdSA(A, phi, mesh, matrices, params, B_mu, p);

% f = @(A) test_Sloc(A, params,matrices,B_mu,phi,p);
% g = @(A) Sloc + Valve_GetdSA(A, phi, mesh, matrices, params, B_mu, p);

for i = 1:2
    if i == 1
        dir = mesh.x;           
    else
        dir = sin(mesh.x + mesh.y);
    end
    err = Diff_Derivatives(f, g, A0, dir);

    fprintf('The relative error = %1.3e\n', err);
end



