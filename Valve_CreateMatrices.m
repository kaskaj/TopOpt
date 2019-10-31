function Valve_CreateMatrices(params, mesh0, refin_level, check)
    
    if nargin < 4
        check = 0;
    end
    
    tol0 = 0.01;   % To refine mesh around the iron
    tol1 = 0.002;  % To refine mesh around the piston
    tol2 = 0.0005; % Size of air where the force on the piston is computed
    
    %% Refine mesh
    
    nodes2coord  = mesh0.nodes2coord;
    bedges2nodes = mesh0.bedges2nodes;
    elems2nodes  = mesh0.elems2nodes;
    tnum         = mesh0.tnum;
    
    for i=1:refin_level
        x_patch = [nodes2coord(elems2nodes(:,1),1), nodes2coord(elems2nodes(:,2),1), nodes2coord(elems2nodes(:,3),1)];
        y_patch = [nodes2coord(elems2nodes(:,1),2), nodes2coord(elems2nodes(:,2),2), nodes2coord(elems2nodes(:,3),2)];
        
        if i <= 3
            % Refine everywhere
            ele_size  = sum(abs(x_patch(:,:) - x_patch(:,[2 3 1])) + abs(y_patch(:,:) - y_patch(:,[2 3 1])),2);
            ii_refine = ele_size >= mean(ele_size);
        elseif i <= 6
            % Refine only around the iron
            ii = min(x_patch,[],2) >= params.x_fe_min - tol0 & max(x_patch,[],2) <= params.x_fe_max + tol0 ...
                & min(y_patch,[],2) >= params.y_fe_min - tol0 & max(y_patch,[],2) <= params.y_fe_max + tol0;
            ele_size  = sum(abs(x_patch(:,:) - x_patch(:,[2 3 1])) + abs(y_patch(:,:) - y_patch(:,[2 3 1])),2);
            ii_refine = ii & ele_size >= mean(ele_size(ii));
        else
            % Refine only around the piston
            x_piston_min = params.x_piston_min;
            x_piston_max = params.x_piston_max;
            y_piston_min = params.y_piston_min;
            y_piston_max = params.y_piston_max;
            
            ii_refine1 = all(abs(x_patch-x_piston_min) <= tol1, 2) & all(y_patch >= y_piston_min-tol1, 2) & all(y_patch <= y_piston_max+tol1, 2);
            ii_refine2 = all(abs(x_patch-x_piston_max) <= tol1, 2) & all(y_patch >= y_piston_min-tol1, 2) & all(y_patch <= y_piston_max+tol1, 2);
            ii_refine3 = all(abs(y_patch-y_piston_min) <= tol1, 2) & all(x_patch >= x_piston_min-tol1, 2) & all(x_patch <= x_piston_max+tol1, 2);
            ii_refine4 = all(abs(y_patch-y_piston_max) <= tol1, 2) & all(x_patch >= x_piston_min-tol1, 2) & all(x_patch <= x_piston_max+tol1, 2);
            ii_refine  = ii_refine1 | ii_refine2 | ii_refine3 | ii_refine4;
        end
        [nodes2coord,bedges2nodes,elems2nodes,tnum] = tridiv2(nodes2coord,bedges2nodes,elems2nodes,tnum,ii_refine);
        
        if i == 3 || i == 6
            [nodes2coord,bedges2nodes,elems2nodes,tnum] = smooth2(nodes2coord,bedges2nodes,elems2nodes,tnum);
        end
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
    
    id_dirichlet = x == max(x) | y == min(y) | y == max(y);
    
    %% Prescribe mu and gamma
    
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
    mu            = zeros(nelement,1);
    mu(mu_nodes)  = params.mu0*params.mur;
    mu(~mu_nodes) = params.mu0;
    mu            = repmat(mu,1,9);
    
    %% Transformation for Triangulation
    
    ii        = zeros(nelement,9);
    jj        = zeros(nelement,9);
    
    edet_aa   = zeros(nelement,1);
    mloc_aa   = zeros(nelement,9);
    slocxx_aa = zeros(nelement,9);
    slocyy_aa = zeros(nelement,9);
    clocx_aa  = zeros(nelement,9);
    clocy_aa  = zeros(nelement,9);
    mloc_aa_plunger  = zeros(nelement,9);
    clocx_aa_plunger = zeros(nelement,9);
    clocy_aa_plunger = zeros(nelement,9);
    
    for k = 1:nelement
        [edet,dFinv] = GenerateTransformation(k,elems2nodes,x,y);
        [~,slocxx,~,slocyy,~,mloc,clocx,clocy] = LocalMatrices(edet,dFinv);
        
        ii(k,:) = repmat(elems2nodes(k,:), 1, 3);
        jj(k,:) = reshape(repmat(elems2nodes(k,:), 3, 1), 1, 9);
        
        edet_aa(k,1)     = edet;
        mloc_aa(k,:)     = mloc(:);
        slocxx_aa(k,:)   = slocxx(:);
        slocyy_aa(k,:)   = slocyy(:);
        clocx_aa(k,:)    = clocx(:);
        clocy_aa(k,:)    = clocy(:);
        
        if x_mid(k) >= params.x_piston_min-tol2 && x_mid(k) <= params.x_piston_max+tol2 && y_mid(k) >= params.y_piston_min-tol2 && y_mid(k) <= params.y_piston_max+tol2
            mloc_aa_plunger(k,:)  = mloc(:);
            clocx_aa_plunger(k,:) = clocx(:);
            clocy_aa_plunger(k,:) = clocy(:);
        end
    end
    
    %% Assemble matrices
    
    Mloc    = sparse(ii(:),jj(:),mloc_aa(:));
    Sloc    = sparse(ii(:),jj(:),(slocxx_aa(:)+slocyy_aa(:)));
    Clocx   = sparse(ii(:),jj(:),clocx_aa(:));
    Clocy   = sparse(ii(:),jj(:),clocy_aa(:));
    Sloc_mu = sparse(ii(:),jj(:),(slocxx_aa(:)+slocyy_aa(:))./mu(:));
    Mloc_plunger  = sparse(ii(:),jj(:),mloc_aa_plunger(:));
    Clocx_plunger = sparse(ii(:),jj(:),clocx_aa_plunger(:));
    Clocy_plunger = sparse(ii(:),jj(:),clocy_aa_plunger(:));
    
    %% Save results
    
    mesh = [];
    mesh.x            = x;
    mesh.y            = y;
    mesh.x_mid        = x_mid;
    mesh.y_mid        = y_mid;
    mesh.npoint       = npoint;
    mesh.nelement     = nelement;
    mesh.nodes2coord  = nodes2coord;
    mesh.bedges2nodes = bedges2nodes;
    mesh.elems2nodes  = elems2nodes;
    mesh.id_dirichlet = id_dirichlet;
    
    file_name_mesh   = fullfile('Valve_Data', sprintf('Mesh%d.mat', refin_level));
    file_name_matrix = fullfile('Valve_Data', sprintf('Matrices%d.mat', refin_level));
    
    save(file_name_mesh,'mesh');
    save(file_name_matrix,'ii','jj','Mloc','Sloc','Clocx','Clocy','Sloc_mu','Mloc_plunger','Clocx_plunger','Clocy_plunger');
    
end

