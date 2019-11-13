function dSA = Valve_GetdSA(A, B_ele, phi, mesh, matrices, params, B_mu, p, mu_fe, dmu_fe)

if nargin < 9
    mu_fe = Valve_GetMu(B_ele,B_mu,params);
    dmu_fe = Valve_GetdMu(B_ele,B_mu,params);
end

mu0 = params.mu0;

A_ele = A(mesh.elems2nodes);

v_y = - 2*(phi.^p).*dmu_fe.*B_ele(:,1) ./ ((1-phi)*mu0 + (phi.^p).*mu_fe).^2; %B_ele(:,1) ->  Clocy
v_x =   2*(phi.^p).*dmu_fe.*B_ele(:,2) ./ ((1-phi)*mu0 + (phi.^p).*mu_fe).^2; %B_ele(:,2) -> -Clocx

hat = [sum(matrices.sloc_aa(:,1:3).*A_ele,2), sum(matrices.sloc_aa(:,4:6).*A_ele,2), sum(matrices.sloc_aa(:,7:9).*A_ele,2)];

hat_v_y = hat.*repmat(v_y,1,3);
hat_v_x = hat.*repmat(v_x,1,3);

clocx_ele_aa = matrices.clocx_ele_aa;
clocy_ele_aa = matrices.clocy_ele_aa;

dsa_y = [hat_v_y(:,1).*clocy_ele_aa(:,1), hat_v_y(:,1).*clocy_ele_aa(:,2), hat_v_y(:,1).*clocy_ele_aa(:,3), ...
    hat_v_y(:,2).*clocy_ele_aa(:,1), hat_v_y(:,2).*clocy_ele_aa(:,2), hat_v_y(:,2).*clocy_ele_aa(:,3), ...
    hat_v_y(:,3).*clocy_ele_aa(:,1), hat_v_y(:,3).*clocy_ele_aa(:,2), hat_v_y(:,3).*clocy_ele_aa(:,3)];
dsa_x = [hat_v_x(:,1).*clocx_ele_aa(:,1), hat_v_x(:,1).*clocx_ele_aa(:,2), hat_v_x(:,1).*clocx_ele_aa(:,3), ...
    hat_v_x(:,2).*clocx_ele_aa(:,1), hat_v_x(:,2).*clocx_ele_aa(:,2), hat_v_x(:,2).*clocx_ele_aa(:,3), ...
    hat_v_x(:,3).*clocx_ele_aa(:,1), hat_v_x(:,3).*clocx_ele_aa(:,2), hat_v_x(:,3).*clocx_ele_aa(:,3)];
    
% dsa_x = dsa_x(:,[1 4 7 2 5 8 3 6 9]);
% dsa_y = dsa_y(:,[1 4 7 2 5 8 3 6 9]);

dSA = sparse(matrices.ii(:),matrices.jj(:),dsa_x(:) + dsa_y(:));

end



