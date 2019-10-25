function CreateMatrices(file_name_mesh,file_name_matrix,nodes2coord,elems2nodes,bedges2nodes,id)

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

for k = 1:nelement
    [edet,dFinv] = GenerateTransformation(k,elems2nodes,x,y);
    [~,slocxx,~,slocyy,~,mloc,clocx,clocy] = LocalMatrices(edet,dFinv);
    
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
end

Mloc  = sparse(ii(:),jj(:),mloc_aa(:));
Sloc  = sparse(ii(:),jj(:),(slocxx_aa(:)+slocyy_aa(:)));
Clocx = sparse(ii(:),jj(:),clocx_aa(:));
Clocy = sparse(ii(:),jj(:),clocy_aa(:));

save(file_name_mesh,'nodes2coord','bedges2nodes','elems2nodes','id','x','y','x_mid','y_mid','npoint','nelement');
save(file_name_matrix,'ii','jj','Mloc','Sloc','Clocx','Clocy','mloc_aa','slocxx_aa','slocyy_aa','clocx_aa','clocy_aa');

end

