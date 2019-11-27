function [A, B] = Motor_GetJ(phi, mesh, matrices, params, model)

id     = ~mesh.id_dirichlet & ~mesh.id_s1 & ~mesh.id_s2 & ~mesh.id_s3;

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

f_j = matrices.Mloc*matrices.J;

%% Compute A

mu_fe   = params.mu0*params.mur;
mu_inv  = 1./((1-phi)*mu_air + (phi.^p).*mu_fe);
mu_inv  = repmat(mu_inv,1,9);
Sloc_mu = sparse(matrices.ii(:),matrices.jj(:),((matrices.sloc_aa(:)).*mu_inv(:)));

% S = [Sloc_mu(id,id),Sloc_mu(id,s1)+Sloc_mu(id,s2);...
%     Sloc_mu(s1,id)+Sloc_mu(s2,id),Sloc_mu(s1,s1)+Sloc_mu(s1,s2)+Sloc_mu(s2,s1)+Sloc_mu(s2,s2)];
% f = [f_j(id);f_j(s1)+f_j(s2)];

% S = [Sloc_mu(id,id),Sloc_mu(id,s1)-Sloc_mu(id,s2);...
%     Sloc_mu(s1,id),Sloc_mu(s1,s1)-Sloc_mu(s1,s2)];
% f = [f_j(id);f_j(s1)+f_j(s2)];
% 
S = [Sloc_mu(id,id),Sloc_mu(id,s1),Sloc_mu(id,s2);...
    Sloc_mu(s1,id)-Sloc_mu(s2,id),Sloc_mu(s1,s1)-Sloc_mu(s2,s1),Sloc_mu(s1,s2)-Sloc_mu(s2,s2);...
    sparse(length(s1),sum(id)), speye(length(s1)), speye(length(s1))];

f = [f_j(id);f_j(s1)+f_j(s2);zeros(length(s1),1)];


% S = [Sloc_mu(id,id),Sloc_mu(id,s1),Sloc_mu(id,s2);...
%     Sloc_mu(s1,id),Sloc_mu(s1,s1),Sloc_mu(s1,s2);...
%     sparse(length(s1),sum(id)), speye(length(s1)), speye(length(s1))];
% S = [Sloc_mu(id,id),Sloc_mu(id,s1),Sloc_mu(id,s2);...
%     Sloc_mu(s2,id),Sloc_mu(s2,s1),Sloc_mu(s2,s2);...
%     sparse(length(s1),sum(id)), speye(length(s1)), speye(length(s1))];
% 
% f = [f_j(id);f_j(s1)+f_j(s2);zeros(length(s1),1)];

qwe = S \ f;

A     = zeros(npoint,1);
A(id) = qwe(1:sum(id));
A(s1) = qwe(sum(id)+1:sum(id)+length(s1));
% A(s2) = qwe(sum(id)+length(s1)+1:end);
A(s2) = -qwe(sum(id)+1:sum(id)+length(s1));

% figure;
% plot(sort(qwe));

%% Compute B

B     = [matrices.Mloc\(matrices.Clocy*A),-matrices.Mloc\(matrices.Clocx*A)];

end

