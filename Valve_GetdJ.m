function [dJdphi,dJdp] = Valve_GetdJ(phi, A, B, B_ele, Sloc_mu, mu_fe, dmu_fe, mesh, matrices, params, model)

id     = ~mesh.id_dirichlet;
npoint = mesh.npoint;

mu_air      = params.mu0;
p           = model.p;

%% Compute beta and gamma

Cp_x  = (-matrices.Clocy_plunger' - matrices.Clocy_plunger)*B(:,1) + (matrices.Clocx_plunger' + matrices.Clocx_plunger)*B(:,2);
Cp_y  = (matrices.Clocx_plunger' + matrices.Clocx_plunger)*B(:,1) + (matrices.Clocy_plunger' + matrices.Clocy_plunger)*B(:,2);
beta  = -(1/params.mu0) * (matrices.Mloc\Cp_x);
gamma = -(1/params.mu0) * (matrices.Mloc\Cp_y);

%% Compute alpha

f         = matrices.Clocy'*beta - matrices.Clocx'*gamma;
alpha     = zeros(npoint,1);

if model.nonlinear == 1
    dSA  = Valve_GetdSA(phi, A, B_ele, mesh, matrices, params, model, mu_fe, dmu_fe);
    dSAf = Sloc_mu(id,id) + dSA(id,id);
    
    alpha(id) = dSAf\f(id);
else
    alpha(id) = Sloc_mu(id,id)\f(id);
end


%% Compute derivative dJdphi

dmu_inv = repmat((mu_air - p.*mu_fe.*phi.^(p-1))./((1-phi)*mu_air + (phi.^p).*mu_fe).^2, 1, 9);
dmu_inv = reshape(dmu_inv', [], 1);

x1 = A(mesh.elems2nodes);
x2 = alpha(mesh.elems2nodes);
y1 = reshape(repmat(x1, 1, 3)', [], 1);
y2 = reshape([repmat(x2(:,1), 1, 3), repmat(x2(:,2), 1, 3), repmat(x2(:,3), 1, 3)]', [], 1);

dJdphi = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* dmu_inv;
dJdphi = sum(reshape(dJdphi, 9, []))';
dJdphi = -dJdphi;

%% Compute derivative dJdp

if nargout > 1
    
    if model.nonlinear == 1
        B_mu	= model.B_mu;
        
        if strcmp(B_mu.approx,'Weibull') == 0
            error('The derivative dJdp can be calculated only for Weibull approximation.');
        end
        
        B0      = (B_ele(:,1).^2) + (B_ele(:,2).^2);
        
        da1 = -(phi.^p).*(B_mu.a_w(3).*exp(-(B0-B_mu.a_w(1)).^B_mu.a_w(2)).*(B_mu.a_w(2).*(B0-B_mu.a_w(1)).^(2*B_mu.a_w(2)-2) - (B_mu.a_w(2)-1).*(B0-B_mu.a_w(1)).^(B_mu.a_w(2)-2)))./((1-phi)*mu_air + (phi.^p).*mu_fe).^2;
        da2 = -(phi.^p).*(B_mu.a_w(3).*exp(-(B0-B_mu.a_w(1)).^B_mu.a_w(2)).*(log(B0-B_mu.a_w(1)).*(B0-B_mu.a_w(1)).^(B_mu.a_w(2)-1) - log(B0-B_mu.a_w(1)).*(B0-B_mu.a_w(1)).^(2*B_mu.a_w(2)-1)))./((1-phi)*mu_air + (phi.^p).*mu_fe).^2;
        da3 = -(phi.^p).*((B0-B_mu.a_w(1)).^(B_mu.a_w(2)-1).*exp(-(B0-B_mu.a_w(1)).^B_mu.a_w(2)))./((1-phi)*mu_air + (phi.^p).*mu_fe).^2;
        
        da1 = reshape(repmat(da1,1,9)', [], 1);
        da2 = reshape(repmat(da2,1,9)', [], 1);
        da3 = reshape(repmat(da3,1,9)', [], 1);
        
        dJdp1 = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* da1;
        dJdp2 = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* da2;
        dJdp3 = reshape(matrices.sloc_aa', [], 1) .* y1 .* y2 .* da3;
        
        dJdp1 = sum(reshape(dJdp1, 9, []))';
        dJdp2 = sum(reshape(dJdp2, 9, []))';
        dJdp3 = sum(reshape(dJdp3, 9, []))';
        
        dJdp = [-sum(dJdp1),-sum(dJdp2),-sum(dJdp3)];
        
    else
        error('The derivative dJdp cannot be calculated for a linear model');        
    end
    
end

end

