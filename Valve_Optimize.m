%% Load data
refin_level = 1;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrices');

%% Prescribe fixed Air and Iron domains
    %tnum: 1 - Air; 2,3 - Coils; 4 - Plunger; 5 - Iron;
    %6 - region under the plunger; 7 - region above the plunger;
    %8 - Channel)
phi = zeros(mesh.nelement,1);

%Air
phi(mesh.tnum == 2) = 0;
phi(mesh.tnum == 3) = 0;
phi(mesh.tnum == 6) = 0;
phi(mesh.tnum == 7) = 0;
phi(mesh.tnum == 8) = 0;

%Iron
phi(mesh.tnum == 4) = 1;

%Initial soloution (Iron)
mu_nodes = (mesh.x_mid >= params.x_init_min) & (mesh.x_mid <= params.x_init_max) & (mesh.y_mid >= params.y_init_min) & (mesh.y_mid <= params.y_init_max);
phi(mu_nodes)  = 1;

