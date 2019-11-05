clear all;

levels_max = 5;

% Set parameters
Valve_SetParameters;

% Load parameters
file_name = fullfile('Valve_Data', 'Param.mat');
load(file_name, 'params');

% Create mesh
Valve_MeshGen(params);

% Load mesh
file_name = fullfile('Valve_Data', 'Mesh0.mat');
load(file_name, 'mesh');

% Create matrices
for levels = 1:levels_max
    fprintf('Creating matrix level %d.\n', levels);
    if levels <= 5
        Valve_CreateMatrices(params, mesh, levels, 1);
    else
        Valve_CreateMatrices(params, mesh, levels, 0);
    end
end

