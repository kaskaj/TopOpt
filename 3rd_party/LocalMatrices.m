function [slocxx,slocyy,mloc,slocx,slocy,slocxy] = LocalMatrices(edet,dFinv)
    dphi   = dFinv'*[-1 1 0; -1 0 1];      
    slocxx = 1/2 * dphi(1,:)'*dphi(1,:) * edet;       % slocxx 
    slocyy = 1/2 * dphi(2,:)'*dphi(2,:) * edet;       % slocyy    
    mloc   = edet*[1 1/2 1/2;1/2 1 1/2;1/2 1/2 1]/12; % mloc        
    slocx  = dphi(1,:)'*edet/2;           % slocx     
    slocy  = dphi(2,:)'*edet/2;           % slocy  
    slocxy = 1/2 * dphi(1,:)'*dphi(2,:) * edet;       % slocxy
end