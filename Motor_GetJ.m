function [T, A, B, B_ele, Sloc_mu, mu_fe, dmu_fe, f] = Motor_GetJ(phi, J, mesh, matrices, params, model, A0)

id     = ~mesh.id_dirichlet & ~mesh.id_s1 & ~mesh.id_s2 & ~mesh.id_s3;
id2    = ~mesh.id_dirichlet & ~mesh.id_s3;

s1 = find(mesh.id_s1);
[~, order] = sort(mesh.y(s1));
s1 = s1(order);

s2 = find(mesh.id_s2);
[~, order] = sort(mesh.x(s2));
s2 = s2(order);

npoint = mesh.npoint;

mu_air = params.mu0;
p      = model.p;

%% Right-hand side

f_j = matrices.Mloc*J;
f = [f_j(id);f_j(s1)+f_j(s2);zeros(length(s1),1)];

%% Compute A

if model.nonlinear == 1
    
    % If initial guess is not provided, compute A
    if nargin < 7 || isempty(A0)
        mu_fe   = params.mu0*params.mur;
        mu_inv  = 1./((1-phi)*mu_air + (phi.^p).*mu_fe);
        mu_inv  = repmat(mu_inv,1,9);
        Sloc_mu = sparse(matrices.ii(:),matrices.jj(:),((matrices.sloc_aa(:)).*mu_inv(:)));
        S       = [Sloc_mu(id,id),Sloc_mu(id,s1),Sloc_mu(id,s2);...
                   Sloc_mu(s1,id)-Sloc_mu(s2,id),Sloc_mu(s1,s1)-Sloc_mu(s2,s1),Sloc_mu(s1,s2)-Sloc_mu(s2,s2);...
                   sparse(length(s1),sum(id)), speye(length(s1)), speye(length(s1))];        

        A_i   = S \ f;

        A     = zeros(npoint,1);
        A(id) =  A_i(1:sum(id));
        A(s1) =  A_i(sum(id)+1:sum(id)+length(s1));
        A(s2) = -A_i(sum(id)+1:sum(id)+length(s1));
    else
        A     = A0;
    end
    B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];
    
    maxsteps = 100;
    
    convergence = zeros(maxsteps,2);    
    
    for i = 1:maxsteps
        
        % Update mu
        [mu_fe, dmu_fe] = Valve_GetMu(B_ele, params, model);

        mu_inv = 1./((1-phi)*params.mu0 + (phi.^p).*mu_fe);
        mu_inv = repmat(mu_inv,1,9);
        
        % Update Sloc
        Sloc_mu = sparse(matrices.ii(:),matrices.jj(:),(matrices.sloc_aa(:)).*mu_inv(:));
        S       = [Sloc_mu(id,id),Sloc_mu(id,s1),Sloc_mu(id,s2);...
                   Sloc_mu(s1,id)-Sloc_mu(s2,id),Sloc_mu(s1,s1)-Sloc_mu(s2,s1),Sloc_mu(s1,s2)-Sloc_mu(s2,s2);...
                   sparse(length(s1),sum(id)), speye(length(s1)), speye(length(s1))];
        
        %Function evaluation and its derivative for S_new and A_old
        SAf       = S*A_i - f;       
          
        dSA  = Valve_GetdSA(phi, A, B_ele, mesh, matrices, params, model, mu_fe, dmu_fe);
        dSAf = S + [dSA(id,id),dSA(id,s1),dSA(id,s2); ...
                        dSA(s1,id)-dSA(s2,id),dSA(s1,s1)-dSA(s2,s1),dSA(s1,s2)-dSA(s2,s2); ...
                        sparse(length(s1),sum(id)), speye(length(s1)), speye(length(s1))];       
        
        % Compute residual
        res = SAf'*matrices.Mloc(id2,id2)*SAf;
        if res <= 1e-8
            break;
        end
        
        % Set step
        if res <= 1e-1
            step = 1;
        else
            step = 0.5;
        end
        
        % Update A and B_ele 
        A_i  =  A_i - step*(dSAf \ SAf);

        A     = zeros(npoint,1);
        A(id) =  A_i(1:sum(id));
        A(s1) =  A_i(sum(id)+1:sum(id)+length(s1));
        A(s2) = -A_i(sum(id)+1:sum(id)+length(s1));
        B_ele  = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];
        
        convergence(i,1) = step;
        convergence(i,2) = mean(A);
        
    end
        fprintf('%d Newton iterations\n',i);
%         hold on;       
%         ax = plotyy(1:(i-1),convergence(1:(i-1),2),1:(i-1),convergence(1:(i-1),1));        
%         title('Convergence plot');
%         xlabel('steps');
%         ylabel(ax(1),'mean(A)');    
%         ylabel(ax(2),'Stepsize');        
    
    if i == maxsteps
        warning('Newton method did not converge.');
    end
else
        
    mu_fe  = params.mu0*params.mur;
    dmu_fe = [];

    mu_inv = 1./((1-phi)*mu_air + (phi.^p).*mu_fe);
    mu_inv = repmat(mu_inv,1,9);
    
    Sloc_mu  = sparse(matrices.ii(:),matrices.jj(:),((matrices.sloc_aa(:)).*mu_inv(:)));    
    S = [Sloc_mu(id,id),Sloc_mu(id,s1),Sloc_mu(id,s2);...
             Sloc_mu(s1,id)-Sloc_mu(s2,id),Sloc_mu(s1,s1)-Sloc_mu(s2,s1),Sloc_mu(s1,s2)-Sloc_mu(s2,s2);...
             sparse(length(s1),sum(id)), speye(length(s1)), speye(length(s1))];        

    A_i = S \ f;
    A     = zeros(npoint,1);
    A(id) =  A_i(1:sum(id));
    A(s1) =  A_i(sum(id)+1:sum(id)+length(s1));
    A(s2) = -A_i(sum(id)+1:sum(id)+length(s1));
end

%% Compute B

B     = [matrices.Mloc\(matrices.Clocy*A),-matrices.Mloc\(matrices.Clocx*A)];
B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];

%% Arrkio's method

r      =  sqrt(mesh.x.^2 + mesh.y.^2);
Br     =  ( B(:,1).*mesh.x + B(:,2).*mesh.y)./r;
Bphi   =  (-B(:,1).*mesh.y + B(:,2).*mesh.x)./r;

ii_nan_r   = find(isnan(Br));
ii_nan_phi = find(isnan(Bphi));

Br(ii_nan_r)     = 0;
Bphi(ii_nan_phi) = 0;

R1_gap0 = params.D3/2;
R2_gap0 = params.D1/2;
d0 = R2_gap0 - R1_gap0;

R1_gap1 = params.D3/2;
R2_gap1 = params.D1/2 - params.d/4;
d1 = R2_gap1 - R1_gap1;

R1_gap2 = params.D3/2;
R2_gap2 = params.D1/2 - params.d/2;
d2 = R2_gap2 - R1_gap2;

R1_gap3 = params.D3/2;
R2_gap3 = params.D1/2 - 3*params.d/4;
d3 = R2_gap3 - R1_gap3;

T0 = params.L*(r'.*Br'*matrices.Mloc_rot0*Bphi)*(1/params.mu0)*(1/d0);
T1 = params.L*(r'.*Br'*matrices.Mloc_rot1*Bphi)*(1/params.mu0)*(1/d1);
T2 = params.L*(r'.*Br'*matrices.Mloc_rot2*Bphi)*(1/params.mu0)*(1/d2);
T3 = params.L*(r'.*Br'*matrices.Mloc_rot3*Bphi)*(1/params.mu0)*(1/d3);

T = [T0,T1,T2,T3];

end

