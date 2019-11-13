clear all;
close all;

%% Load data
refin_level = 4;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, 'B_mu'),'B_mu');
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

id     = ~mesh.id_dirichlet;

%% Newton–Raphson

p = 1;
coil = 1;
nonlinear = 1;
steps = 100;

%Initial guess
mu_fe = params.mu0*params.mur*ones(mesh.nelement,1);
[~, A, ~, B_ele, Sloc, f] = Valve_GetJ(phi, mesh, matrices, params, p, coil, nonlinear, mu_fe);

SAf = zeros(mesh.npoint,1);

for i = 1:steps
    
    %New mu and dmu
    mu_fe = Valve_GetMu(B_ele,B_mu,params);
    mu_inv = 1./((1-phi)*params.mu0 + (phi.^p).*mu_fe);
    mu_inv = repmat(mu_inv,1,9);    
    dmu_fe = Valve_GetdMu(B_ele,B_mu,params);
    
    %New Sloc
    Sloc  = sparse(matrices.ii(:),matrices.jj(:),(matrices.sloc_aa(:)).*mu_inv(:));
    
    %Function evaluation and its derivative for S_new and A_old
    
    SAf(id)  = Sloc(id,id)*A(id) - f(id);
    SAf(~id) = 0;
    
    dSA = Valve_GetdSA(A, B_ele, phi, mesh, matrices, params, B_mu, p, mu_fe, dmu_fe);
    dSAf = Sloc(id,id) + dSA(id,id);
    
    %New A, new B_ele  
    
    res = SAf'*matrices.Mloc*SAf;
    
    if res <= 1e-8
        break;
    end
    
    if res <= 1e-1
        step = 1;
    else
        step = 0.5;
    end

    A(id)  =  A(id) - step*(dSAf \ SAf(id));
    A(~id) = 0;
    B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];

end

[F, A, B, ~, ~, ~] = Valve_GetJ(phi, mesh, matrices, params, p, coil, nonlinear, mu_fe);

fprintf('Force for nonlinear model: Fy = %d\n',F);

%Plot
% ele = delaunay(mesh.x_mid,mesh.y_mid);
% PlotData(mesh.x_mid,mesh.y_mid,ele,B_ele(:,1))
% Valve_PlotEdges(params,max(B_ele(:,1)))


