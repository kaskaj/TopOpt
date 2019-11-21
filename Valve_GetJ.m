function [F, A, B, B_ele, Sloc_mu, mu_fe, dmu_fe, f] = Valve_GetJ(phi, mesh, matrices, params, model, A0)

id     = ~mesh.id_dirichlet;
npoint = mesh.npoint;

mu_air = params.mu0;
p      = model.p;

%% Right-hand side

if model.coil == 1    %Turn on the current
    f = matrices.Mloc*matrices.J + (1/params.mu0)*matrices.Clocy*matrices.Br;
else            %Turn off the current
    f = (1/params.mu0)*matrices.Clocy*matrices.Br;
end

%% Compute A

if model.nonlinear == 1
    
    % If initial guess is not provided, compute A
    if nargin < 6 || isempty(A0)
        mu_fe   = params.mu0*params.mur;        
        mu_inv  = 1./((1-phi)*mu_air + (phi.^p).*mu_fe);
        mu_inv  = repmat(mu_inv,1,9);
        Sloc_mu = sparse(matrices.ii(:),matrices.jj(:),((matrices.sloc_aa(:)).*mu_inv(:)));
        
        A     = zeros(npoint,1);
        A(id) = Sloc_mu(id,id) \ f(id);
    else
        A     = A0;
    end
    B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];
    
    SAf = zeros(mesh.npoint,1);
    
    maxsteps = 100;
    for i = 1:maxsteps
        
        % Update mu
        [mu_fe, dmu_fe] = Valve_GetMu(B_ele, params, model);

        mu_inv = 1./((1-phi)*params.mu0 + (phi.^p).*mu_fe);
        mu_inv = repmat(mu_inv,1,9);
        
        % Update Sloc
        Sloc_mu = sparse(matrices.ii(:),matrices.jj(:),(matrices.sloc_aa(:)).*mu_inv(:));
        
        %Function evaluation and its derivative for S_new and A_old
        SAf(id)  = Sloc_mu(id,id)*A(id) - f(id);
        SAf(~id) = 0;
        
        dSA  = Valve_GetdSA(phi, A, B_ele, mesh, matrices, params, model, mu_fe, dmu_fe);
        dSAf = Sloc_mu(id,id) + dSA(id,id);
        
        % Compute residual
        res = SAf'*matrices.Mloc*SAf;
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
        A(id)  =  A(id) - step*(dSAf \ SAf(id));
        A(~id) = 0;
        B_ele  = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];
        
    end
        %fprintf('%d Newton iterations\n',i);
    
    if i == maxsteps
        warning('Newton method did not converge.');
    end
else
    mu_fe  = params.mu0*params.mur;
    dmu_fe = [];

    mu_inv = 1./((1-phi)*mu_air + (phi.^p).*mu_fe);
    mu_inv = repmat(mu_inv,1,9);
    
    Sloc_mu  = sparse(matrices.ii(:),matrices.jj(:),((matrices.sloc_aa(:)).*mu_inv(:)));
    
    A     = zeros(npoint,1);    
    A(id) = Sloc_mu(id,id) \ f(id);

end

%% Compute B

B     = [matrices.Mloc\(matrices.Clocy*A),-matrices.Mloc\(matrices.Clocx*A)];
B_ele = [matrices.Clocy_ele*A,-matrices.Clocx_ele*A];

%% Compute F

F_x_aux = B(:,1)'*matrices.Clocx_plunger*B(:,1) - B(:,2)'*matrices.Clocx_plunger*B(:,2) + B(:,2)'*matrices.Clocy_plunger*B(:,1) + B(:,1)'*matrices.Clocy_plunger*B(:,2);
F_y_aux = -B(:,1)'*matrices.Clocy_plunger*B(:,1) + B(:,2)'*matrices.Clocy_plunger*B(:,2) + B(:,2)'*matrices.Clocx_plunger*B(:,1) + B(:,1)'*matrices.Clocx_plunger*B(:,2);
F       = 1/params.mu0*[F_x_aux; F_y_aux];

F = F(2);   %force in y direction
F = -F;

% F_mod = zeros(4,2);
% for i=1:4
%     if i==1
%         ii = y == params.y_piston_min & x >= params.x_piston_min & x <= params.x_piston_max;
%     elseif i==2
%         ii = x == params.x_piston_max & y >= params.y_piston_min & y <= params.y_piston_max;
%     elseif i==3
%         ii = y == params.y_piston_max & x >= params.x_piston_min & x <= params.x_piston_max;
%     else
%         ii = x == params.x_piston_min & y >= params.y_piston_min & y <= params.y_piston_max;
%     end
%     x_p  = x(ii);
%     y_p  = y(ii);
%     Bx_p = B(ii,1);
%     By_p = B(ii,2);
%
%     if i==1 || i==3
%         [~,jj] = sort(x_p);
%     else
%         [~,jj] = sort(y_p);
%     end
%     x_p    = x_p(jj);
%     y_p    = y_p(jj);
%     Bx_p   = Bx_p(jj);
%     By_p   = By_p(jj);
%
%     T_p  = {1/2*(Bx_p.*Bx_p-By_p.*By_p), Bx_p.*By_p; Bx_p.*By_p, 1/2*(By_p.*By_p-Bx_p.*Bx_p)};
%     if i==1
%         T_pn = [-T_p{1,2}, -T_p{2,2}];
%     elseif i==2
%         T_pn = [T_p{1,1}, T_p{2,1}];
%     elseif i==3
%         T_pn = [T_p{1,2}, T_p{2,2}];
%     else
%         T_pn = [-T_p{1,1}, -T_p{2,1}];
%     end
%     xy_diff = diff([x_p, y_p]);
%     xy_norm = sqrt(xy_diff(:,1).^2 + xy_diff(:,2).^2);
%
%     T_mean  = 0.5*(T_pn(1:end-1,:) + T_pn(2:end,:));
%
%     F_mod(i,:) = sum(repmat(xy_norm,1,2).*T_mean)';
% end
% F_mod = F_mod / params.mu0;
% F_mod = sum(F_mod)';

end

