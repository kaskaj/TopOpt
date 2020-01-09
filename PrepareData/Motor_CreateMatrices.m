function Motor_CreateMatrices(params, mesh0, refin_level, check)

refin_level = sort(refin_level);

if nargin < 4
    check = 0;
end

nodes2coord  = mesh0.nodes2coord;
bedges2nodes = mesh0.bedges2nodes;
elems2nodes  = mesh0.elems2nodes;
tnum         = mesh0.tnum;

for i_level=1:max(refin_level)
    
     fprintf('Creating mesh level %d.\n', i_level);
        
        %% Refine mesh
                  
        x = nodes2coord(:,1);
        y = nodes2coord(:,2);
        
        x_mid = (1/3)*(nodes2coord(elems2nodes(:,1),1) + nodes2coord(elems2nodes(:,2),1) + nodes2coord(elems2nodes(:,3),1));
        y_mid = (1/3)*(nodes2coord(elems2nodes(:,1),2) + nodes2coord(elems2nodes(:,2),2) + nodes2coord(elems2nodes(:,3),2));
        
               
        if i_level > 2
            tol4 = 2e-3;
            ii_refine = x_mid >= sqrt((params.D3/2 - tol4)^2 - y_mid.^2) & ...
                        x_mid <= sqrt((params.D1/2 + tol4)^2 - y_mid.^2);            
        else
            ii_refine = ismember(tnum,3)| ismember(tnum,4) | ismember(tnum,2);
        end
     
        [nodes2coord,bedges2nodes,elems2nodes,tnum] = tridiv2(nodes2coord,bedges2nodes,elems2nodes,tnum,ii_refine);
        [nodes2coord,bedges2nodes,elems2nodes,tnum] = smooth2(nodes2coord,bedges2nodes,elems2nodes,tnum);
    
    
    if any(i_level == refin_level)
        
        fprintf('Creating matrix level %d.\n', i_level);
        
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
        
        tol1 = 1e-5;
        tol2 = 1e-10;
        tol3 = 1e-7;            
        R1 = params.D2/2; 
        R2 = params.D3/2; 
        
        id_dirichlet = (y >= sqrt(R1^2 - tol1 - x.^2))  | (x - tol1 >= sqrt(R1^2 - y.^2));              
        id_s1 = (x <= min(x)+tol2);
        id_s2 = (y <= min(y)+tol2);  
        id_cloc_rot  = x_mid <= sqrt((R2+tol3)^2 - y_mid.^2);  
        
        %Deselect point 0
        id_s3 = x == 0 & y==0;        
        id_s1(id_s3) = 0;
        id_s2(id_s3) = 0;
        
        %Deselect common points
        id_s1(id_dirichlet) = 0;
        id_s2(id_dirichlet) = 0;
        
        %Plot boundaries
        figure;
        plot(x,y,'o');
        hold on;
        plot(x(id_dirichlet),y(id_dirichlet),'ro');
        
        figure;
        plot(x,y,'o');
        hold on;
        plot(x(id_s1),y(id_s1),'ro');
        
        figure;
        plot(x,y,'o');
        hold on;
        plot(x(id_s2),y(id_s2),'ro');           
        
        figure;
        plot(x_mid,y_mid,'o');
        hold on;
        plot(x_mid(id_cloc_rot),y_mid(id_cloc_rot),'ro');
        Motor_PlotEdges(params,1);
                
        %% Transformation for Triangulation
        
        ii        = zeros(nelement,9);
        jj        = zeros(nelement,9);
        
        ii_ele    = zeros(nelement,3);
        jj_ele    = zeros(nelement,3);
        
        edet_aa   = zeros(nelement,1);
        mloc_aa   = zeros(nelement,9);
        slocxx_aa = zeros(nelement,9);
        slocyy_aa = zeros(nelement,9);
        clocx_aa  = zeros(nelement,9);
        clocy_aa  = zeros(nelement,9);
        clocx_ele_aa  = zeros(nelement,3);
        clocy_ele_aa  = zeros(nelement,3);
        clocx_aa_rotor = zeros(nelement,9);
        clocy_aa_rotor = zeros(nelement,9);
        
        for k = 1:nelement
            [edet,dFinv,Cinv] = GenerateTransformation(k,elems2nodes,x,y);
            [~,slocxx,~,slocyy,~,mloc,clocx,clocy,clocx_ele,clocy_ele] = LocalMatrices(edet,dFinv,Cinv);
            
            ii(k,:) = repmat(elems2nodes(k,:), 1, 3);
            jj(k,:) = reshape(repmat(elems2nodes(k,:), 3, 1), 1, 9);
            
            ii_ele(k,:) = [k,k,k];
            jj_ele(k,:) = elems2nodes(k,:);
            
            edet_aa(k,1)     = edet;
            mloc_aa(k,:)     = mloc(:);
            slocxx_aa(k,:)   = slocxx(:);
            slocyy_aa(k,:)   = slocyy(:);
            clocx_aa(k,:)    = clocx(:);
            clocy_aa(k,:)    = clocy(:);
            clocx_ele_aa(k,:)    = clocx_ele(:);
            clocy_ele_aa(k,:)    = clocy_ele(:);                                 
            
            if  x_mid(k) <= sqrt(R2^2 + tol3 - y_mid(k)^2)                 
                clocx_aa_rotor(k,:) = clocx(:);
                clocy_aa_rotor(k,:) = clocy(:);
            end  
            
        end
        %% Assemble matrices
        
        Mloc    = sparse(ii(:),jj(:),mloc_aa(:));
        Sloc    = sparse(ii(:),jj(:),(slocxx_aa(:)+slocyy_aa(:)));
        Clocx   = sparse(ii(:),jj(:),clocx_aa(:));
        Clocy   = sparse(ii(:),jj(:),clocy_aa(:));
        Clocx_ele = sparse(ii_ele(:),jj_ele(:),clocx_ele_aa(:));
        Clocy_ele = sparse(ii_ele(:),jj_ele(:),clocy_ele_aa(:));
        Clocx_rotor = sparse(ii(:),jj(:),clocx_aa_rotor(:));
        Clocy_rotor = sparse(ii(:),jj(:),clocy_aa_rotor(:));
        
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
        mesh.tnum         = tnum;
        mesh.id_dirichlet = id_dirichlet;
        mesh.id_s1        = id_s1;
        mesh.id_s2        = id_s2;
        mesh.id_s3        = id_s3;
        
        matrices = [];
        matrices.ii            = ii;
        matrices.jj            = jj;
        matrices.Mloc          = Mloc;
        matrices.mloc_aa       = mloc_aa;
        matrices.Sloc          = Sloc;
        matrices.sloc_aa       = slocxx_aa + slocyy_aa;
        matrices.Clocx         = Clocx;
        matrices.Clocx_ele     = Clocx_ele;
        matrices.clocx_ele_aa  = clocx_ele_aa;
        matrices.Clocy         = Clocy;
        matrices.Clocy_ele     = Clocy_ele;
        matrices.clocy_ele_aa  = clocy_ele_aa;
        matrices.Clocx_rotor   = Clocx_rotor;
        matrices.Clocy_rotor   = Clocy_rotor;

        
        file_name_mesh   = fullfile('Motor_Data', sprintf('Mesh%d.mat', i_level));
        file_name_matrix = fullfile('Motor_Data', sprintf('Matrices%d.mat', i_level));
        
        save(file_name_mesh,'mesh');
        save(file_name_matrix,'matrices');
        
    end
    
end

end

