function SAid = test_Sloc_id(A, params, mesh, matrices, B_mu, phi, p, Sloc0)

id     = ~mesh.id_dirichlet;
S0A = zeros(mesh.npoint,1);

SA = test_Sloc(A, params, mesh, matrices, B_mu, phi, p);

S0A(id)= Sloc0(id,id)*A(id);

SAid = SA - S0A;

end

