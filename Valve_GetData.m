function [A,B] = Valve_GetData(mesh,matrix,params,Sloc_mu,plot)

if nargin < 5
    plot = 0;
end

%% Extract parameters

id          = ~mesh.id_dirichlet;
npoint      = mesh.npoint;
elems2nodes = mesh.elems2nodes;

x = mesh.nodes2coord(:,1);
y = mesh.nodes2coord(:,2);

%% Solve the system

f     = matrix.Mloc*matrix.J;
A     = zeros(npoint,1);
A(id) = Sloc_mu(id,id) \ f(id);
B     = [matrix.Mloc\(matrix.Clocy*A),-matrix.Mloc\(matrix.Clocx*A)];

%% Plot fields
if plot
    field = A;
    PlotData(x, y, elems2nodes, field);
    Valve_PlotEdges(params, max(field));
    
    field = B(:,1);
    PlotData(x, y, elems2nodes, field);
    Valve_PlotEdges(params, max(field));
end

end




