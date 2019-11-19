function [mu_fe, dmu_fe] = Valve_Newton(A, B_ele, phi, p, f, mesh, matrices, params, B_mu, maxsteps)

SAf = zeros(mesh.npoint,1);
id     = ~mesh.id_dirichlet;

for i = 1:maxsteps
    
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

if i == maxsteps
    warning('Newton did not converge.');
end

i

end

