clear all;

levels_max = 4;

% Set parameters
Valve_SetParameters;

% Load parameters
file_name = fullfile('../Valve_Data', 'Param.mat');
load(file_name, 'params');

% Create mesh
Valve_MeshGen(params);

% Load mesh
file_name = fullfile('../Valve_Data', 'Mesh0.mat');
load(file_name, 'mesh');

% Create matrices
Valve_CreateMatrices(params, mesh, 1:levels_max, 1);
