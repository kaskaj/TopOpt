load(fullfile('Valve_Data', 'Mesh7'));
load(fullfile('Valve_Data', 'Param'))

filename = 'Valve_Data/ValveNeumann.mph';
model = mphopen(filename);
meshstats = mphxmeshinfo(model);

elems2nodes1 = meshstats.elements.tri.dofs' + 1; %Comsol is indexing from 0
nodes2coord1 = meshstats.dofs.coords';

elems2nodes2 = mesh.elems2nodes;
nodes2coord2 = mesh.nodes2coord;

tol = 0.003;
zoom = [params.x_piston_min - tol, params.x_piston_max + tol; params.y_piston_min - tol, params.y_piston_max + tol];

filename = 'Valve_Results/Mesh.png';
meshname = 'Mesh';
SubPlotMesh(elems2nodes1,nodes2coord1,elems2nodes2,nodes2coord2,params,zoom,filename,meshname,1)
