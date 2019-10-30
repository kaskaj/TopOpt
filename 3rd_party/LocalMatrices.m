function [slocx,slocxx,slocy,slocyy,slocxy,mloc,clocx,clocy] = LocalMatrices(edet,dFinv)
    
    dphi   = dFinv'*[-1 1 0; -1 0 1];
    
    slocxx = 1/2 * dphi(1,:)'*dphi(1,:) * edet;       % slocxx
    slocx  = dphi(1,:)*edet/2;                        % slocx
    slocyy = 1/2 * dphi(2,:)'*dphi(2,:) * edet;       % slocyy
    slocy  = dphi(2,:)*edet/2;                        % slocy
    slocxy = 1/2 * dphi(1,:)'*dphi(2,:) * edet;       % slocxy
    mloc   = edet*[1 1/2 1/2;1/2 1 1/2;1/2 1/2 1]/12; % mloc
    clocx  = [slocx;slocx;slocx]/3;                   % clocx
    clocy  = [slocy;slocy;slocy]/3;                   % clocy
end

