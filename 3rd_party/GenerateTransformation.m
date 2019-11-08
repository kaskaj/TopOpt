function [edet,dFinv,Cinv] = GenerateTransformation(k,e2p,x,y)    
    dx1   = x(e2p(k,2))-x(e2p(k,1));
    dy1   = y(e2p(k,2))-y(e2p(k,1));
    dx2   = x(e2p(k,3))-x(e2p(k,1));
    dy2   = y(e2p(k,3))-y(e2p(k,1));
    edet  = dx1.*dy2 - dx2.*dy1;           % Determinant of the Jacobian
    C     = [1, x(e2p(k,1)), y(e2p(k,1));
             1, x(e2p(k,2)), y(e2p(k,2));
             1, x(e2p(k,3)), y(e2p(k,3));
             ];
    Cinv = inv(C);      
    dFinv = 1/edet*[dy2, -dx2; -dy1, dx1]; % Inverse Jacobian
end

