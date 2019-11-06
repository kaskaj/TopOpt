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

%% Select test mu

% Select whole Fe region + air
mu_nodes_select    = (x_mid >= params.x_fe_min) & (x_mid <= params.x_fe_max) & (y_mid >= params.y_fe_min) & (y_mid <= params.y_fe_max);

% Deselect conductors and gaps
mu_nodes_deselect1 = (x_mid >= params.x_c_min) & (x_mid <= params.x_c_max) & (y_mid >= params.y_c_min) & (y_mid <= params.y_c_max);
mu_nodes_deselect2 = (x_mid >= params.x_c_min) & (x_mid <= params.x_c_max) & (y_mid >= -params.y_c_max) & (y_mid <= -params.y_c_min);
mu_nodes_deselect3 = (x_mid >= params.x_gap_min) & (x_mid <= params.x_gap_max) & (y_mid >= params.y_gap_min) & (y_mid <= params.y_gap_max);
mu_nodes_deselect4 = (x_mid >= params.x_t1_min) & (x_mid <= params.x_t1_max) & (y_mid >= params.y_t1_min) & (y_mid <= params.y_t1_max);
mu_nodes_deselect5 = (x_mid >= params.x_t2_min) & (x_mid <= params.x_t2_max) & (y_mid >= params.y_t2_min) & (y_mid <= params.y_t2_max);

% Combine them
mu_nodes           = mu_nodes_select & ~mu_nodes_deselect1 & ~mu_nodes_deselect2 & ~mu_nodes_deselect3 & ~mu_nodes_deselect4 & ~mu_nodes_deselect5;

% Insert values
phi            = zeros(mesh.nelement,1);
phi(mu_nodes)  = 1;

%% Compute derivatives
p = 2;

[F, A, B, Sloc_mu] = Valve_GetJ(phi, mesh, matrices, params, p);
dJ = Valve_GetdJ(phi, Sloc_mu, A, B, mesh, matrices, params, p);

f = @(phi) Valve_GetJ(phi, mesh, matrices, params, p);
g = @(x) Valve_GetdJ(phi, Sloc_mu, A, B, mesh, matrices, params, p);

%Direction

%dir = dJ;
dir = sin(mesh.x_mid + mesh.y_mid); 

test = Diff_Derivatives(f, g, phi, dir);

test
