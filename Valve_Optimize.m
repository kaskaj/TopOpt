clear all;

%% Load Data
refin_level = 1;

folder_name = 'Valve_Data';
load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');
load(fullfile(folder_name, sprintf('Matrices%d.mat', refin_level)), 'matrix');

nelement = size(mesh.elems2nodes,1);

mu1 = params.mu0;
mu2 = params.mu0*params.mur;

%% phi -> mu, dmu

phi = ones(nelement,1);
%phi = rand(nelement,1);

mu_inv = 1./((1-phi)*mu1 + phi*mu2);
mu_inv = repmat(mu_inv,1,9);

dmu_inv = (mu1 - mu2)./((1-phi)*mu1 + phi*mu2).^2;
dmu_inv = repmat(dmu_inv,1,9);

%% New Sloc_mu and dSloc_mu

Sloc_mu  = sparse(matrix.ii(:),matrix.jj(:),(matrix.sloc_aa(:)).*mu_inv(:));
dSloc_mu = sparse(matrix.ii(:),matrix.jj(:),(matrix.sloc_aa(:)).*dmu_inv(:));

%% Soloution for Sloc_mu
plot = 1;
[A,B] = Valve_GetData_test(mesh,matrix,params,Sloc_mu,plot);

%% Compute derivative

Cp_x  = (-matrix.Clocy_plunger' - matrix.Clocy_plunger)*B(:,1) + (matrix.Clocx_plunger' + matrix.Clocx_plunger)*B(:,2);
Cp_y  = (matrix.Clocx_plunger' + matrix.Clocx_plunger)*B(:,1) + (matrix.Clocx_plunger' + matrix.Clocx_plunger)*B(:,2);
beta  = -(1/params.mu0) * (matrix.Mloc_plunger\Cp_x);
gamma = -(1/params.mu0) * (matrix.Mloc_plunger\Cp_y);
alpha = Sloc_mu\(matrix.Clocy'*beta - matrix.Clocx'*gamma);

dJ = alpha'*dSloc_mu*A;

dJ


