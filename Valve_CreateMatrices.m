function Valve_CreateMatrices(levels,check)

load(fullfile('Valve_Data', 'param'));
load(fullfile('Valve_Data', 'mesh0'));

if nargin < 2
    check = 0;
end       

%% Refine mesh

for i=1:levels
    [nodes2coord,elems2nodes,bedges2nodes] = refinement_uniform_2D(nodes2coord,elems2nodes,bedges2nodes);
end

%% Get data

x = nodes2coord(:,1);
y = nodes2coord(:,2);

x_mid = (1/3)*(nodes2coord(elems2nodes(:,1),1) + nodes2coord(elems2nodes(:,2),1) + nodes2coord(elems2nodes(:,3),1));
y_mid = (1/3)*(nodes2coord(elems2nodes(:,1),2) + nodes2coord(elems2nodes(:,2),2) + nodes2coord(elems2nodes(:,3),2));

nelement = size(elems2nodes,1);
npoint   = length(x);

if check    
    for k=1:nelement
        xy = nodes2coord(elems2nodes(k,:),:);
        score = (xy(3,1)-xy(2,1))*(xy(3,2)+xy(2,2)) + (xy(2,1)-xy(1,1))*(xy(2,2)+xy(1,2)) + (xy(1,1)-xy(3,1))*(xy(1,2)+xy(3,2));
        if score >= 0
            error('Some element is represented clockwise');
        end
    end    
end

%% Set boundary conditions

idp = zeros(size(x));
idp(x == min(x) | x == max(x)) = 1;
idp(y == min(y) | y == max(y)) = 1;

id     = ~idp(1:npoint) == 1;

%% Prescribe mu and gamma

mu    = zeros(nelement,1);
gamma = zeros(nelement,1);

%Select whole Fe region + air

mu_nodes_select    = (x_mid >= x_fe_min) & (x_mid <= x_fe_max) & (y_mid >= y_fe_min) & (y_mid <= y_fe_max);

%Deselect conductors and gaps

mu_nodes_deselect1 = (x_mid >= x_c_min) & (x_mid <= x_c_max) & (y_mid >= y_c_min) & (y_mid <= y_c_max);
mu_nodes_deselect2 = (x_mid >= x_c_min) & (x_mid <= x_c_max) & (y_mid >= -y_c_max) & (y_mid <= -y_c_min);
mu_nodes_deselect3 = (x_mid >= x_gap_min) & (x_mid <= x_gap_max) & (y_mid >= y_gap_min) & (y_mid <= y_gap_max);
mu_nodes_deselect4 = (x_mid >= x_t_min) & (x_mid <= x_t_max) & (y_mid >= y_t_min) & (y_mid <= y_t_max);
mu_nodes_deselect5 = (x_mid >= x_t2_min) & (x_mid <= x_t2_max) & (y_mid >= y_t2_min) & (y_mid <= y_t2_max);

mu_nodes = mu_nodes_select & ~mu_nodes_deselect1 & ~mu_nodes_deselect2 & ~mu_nodes_deselect3 & ~mu_nodes_deselect4 & ~mu_nodes_deselect5;

mu(mu_nodes)  = mu0*mur;
mu(~mu_nodes) = mu0;

mu    = repmat(mu,1,9);

file_name_mesh   = fullfile('Valve_Data', sprintf('mesh%d.mat', levels));
file_name_matrix = fullfile('Valve_Data', sprintf('Matrices%d.mat', levels));

CreateMatrices(file_name_mesh,file_name_matrix,nodes2coord,elems2nodes,bedges2nodes,id)

load(file_name_matrix, 'ii', 'jj', 'slocxx_aa', 'slocyy_aa', 'mloc_aa');

Sloc_mu    = sparse(ii(:),jj(:),(slocxx_aa(:)+slocyy_aa(:))./mu(:));

save(file_name_matrix, 'Sloc_mu', '-append');

end

