function CreateMatrices(file_name_mesh,file_name_matrix,nodes2coord,elems2nodes,bedges2nodes,id,mu)

    load(fullfile('Valve_Data', 'param'));
    
x = nodes2coord(:,1);
y = nodes2coord(:,2);

x_mid = (1/3)*(nodes2coord(elems2nodes(:,1),1) + nodes2coord(elems2nodes(:,2),1) + nodes2coord(elems2nodes(:,3),1));
y_mid = (1/3)*(nodes2coord(elems2nodes(:,1),2) + nodes2coord(elems2nodes(:,2),2) + nodes2coord(elems2nodes(:,3),2));

nelement = size(elems2nodes,1);
npoint   = length(x);

%% Transformation for Triangulation

ii        = zeros(nelement,9);
jj        = zeros(nelement,9);

edet_aa   = zeros(nelement,1);
mloc_aa   = zeros(nelement,9);
slocxx_aa = zeros(nelement,9);
slocyy_aa = zeros(nelement,9);
clocx_aa  = zeros(nelement,9);
clocy_aa  = zeros(nelement,9);
glocx_aa  = zeros(nelement,9);
glocy_aa  = zeros(nelement,9);
mloc_aa_plunger   = zeros(nelement,9);
clocx_aa_plunger = zeros(nelement,9);
clocy_aa_plunger = zeros(nelement,9);
glocx_aa_plunger = zeros(nelement,9);
glocy_aa_plunger = zeros(nelement,9);

for k = 1:nelement
    [edet,dFinv] = GenerateTransformation(k,elems2nodes,x,y);
    [~,slocxx,~,slocyy,slocxy,mloc,clocx,clocy,glocx,glocy] = LocalMatrices(edet,dFinv);
    
    e2pRow1 = repmat(elems2nodes(k,:), 1, 3);
    e2pRow2 = kron(elems2nodes(k,:), ones(1, 3));
    ii(k,:) = e2pRow1;
    jj(k,:) = e2pRow2;
    
    edet_aa(k,1)     = edet;
    mloc_aa(k,:)     = mloc(:);
    slocxx_aa(k,:)   = slocxx(:);
    slocyy_aa(k,:)   = slocyy(:);
    clocx_aa(k,:)    = clocx(:);
    clocy_aa(k,:)    = clocy(:);
    glocx_aa(k,:)    = glocx(:);
    glocy_aa(k,:)    = glocy(:);

    if x_mid(k) >= x_channel && x_mid(k) <= x_channel + x1 && y_mid(k) >= -y3/2-2*t-y2+y4+t && y_mid(k) <= -y3/2-2*t-y2+y4+t+y_piston
        mloc_aa_plunger(k,:)   = mloc(:);
        clocx_aa_plunger(k,:) = clocx(:);
        clocy_aa_plunger(k,:) = clocy(:);
        glocx_aa_plunger(k,:)    = glocx(:);
        glocy_aa_plunger(k,:)    = glocy(:);
    end
end

Mloc  = sparse(ii(:),jj(:),mloc_aa(:));
Sloc  = sparse(ii(:),jj(:),(slocxx_aa(:)+slocyy_aa(:)));
Clocx = sparse(ii(:),jj(:),clocx_aa(:));
Clocy = sparse(ii(:),jj(:),clocy_aa(:));
Glocx = sparse(ii(:),jj(:),glocx_aa(:));
Glocy = sparse(ii(:),jj(:),glocy_aa(:));
Sloc_mu = sparse(ii(:),jj(:),(slocxx_aa(:)+slocyy_aa(:))./mu(:));
Mloc_plunger   = sparse(ii(:),jj(:),mloc_aa_plunger(:));
Clocx_plunger = sparse(ii(:),jj(:),clocx_aa_plunger(:));
Clocy_plunger = sparse(ii(:),jj(:),clocy_aa_plunger(:));
Glocx_plunger = sparse(ii(:),jj(:),glocx_aa_plunger(:));
Glocy_plunger = sparse(ii(:),jj(:),glocy_aa_plunger(:));
% Slocyy_plunger = sparse(ii(:),jj(:),slocyy_aa_plunger(:));

save(file_name_mesh,'nodes2coord','bedges2nodes','elems2nodes','id','x','y','x_mid','y_mid','npoint','nelement');
save(file_name_matrix,'ii','jj','Mloc','Sloc','Clocx','Clocy','Glocx','Glocy','Sloc_mu','Mloc_plunger','Clocx_plunger','Clocy_plunger','Glocx_plunger','Glocy_plunger');

end

