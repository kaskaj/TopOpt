function [F, A, B, Sloc_mu] = Valve_GetJ(phi, mesh, matrices, params, p, coil)
    
    if nargin < 6
        coil = 1;
    end
    
    id     = ~mesh.id_dirichlet;
    npoint = mesh.npoint;
    
    mu1 = params.mu0;
    mu2 = params.mu0*params.mur;
    
    %% phi -> mu, dmu
    
    mu_inv = 1./((1-phi)*mu1 + (phi.^p)*mu2);
    mu_inv = repmat(mu_inv,1,9);
    
    Sloc_mu  = sparse(matrices.ii(:),matrices.jj(:),(matrices.sloc_aa(:)).*mu_inv(:));
    
    %% Solve the system
    if coil == 1    %Turn on the current
        f = matrices.Mloc*matrices.J + (1/params.mu0)*matrices.Clocy*matrices.Br;
    else            %Turn off the current
        f = (1/params.mu0)*matrices.Clocy*matrices.Br;
    end
    A     = zeros(npoint,1);
    A(id) = Sloc_mu(id,id) \ f(id);
    B     = [matrices.Mloc\(matrices.Clocy*A),-matrices.Mloc\(matrices.Clocx*A)];
    
    
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

