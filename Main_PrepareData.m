clear all;

addpath(genpath('3rd_party'));
addpath('PrepareData');

levels_max = 3;

% Set parameters
Motor_SetParameters;

% Load parameters
% file_name = fullfile('Valve_Data', 'Param.mat');
file_name = fullfile('Motor_Data', 'Param.mat');
load(file_name, 'params');

% Create mesh
% Valve_MeshGen(params);
Motor_MeshGen(params);

% Load mesh
% file_name = fullfile('Valve_Data', 'Mesh0.mat');
file_name = fullfile('Motor_Data', 'Mesh0.mat');
load(file_name, 'mesh');

% Create matrices
% Valve_CreateMatrices(params, mesh, 1:levels_max, 1);
Motor_CreateMatrices(params, mesh, 1:levels_max, 1);

