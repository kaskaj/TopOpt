refin_level = 4;
save = 1;
tol = 0.001;

folder_name = 'Valve_Data';

load(fullfile(folder_name, 'Param'), 'params');
load(fullfile(folder_name, sprintf('Mesh%d.mat', refin_level)), 'mesh');

x = mesh.x;
y = mesh.y;
elems2nodes = mesh.elems2nodes;

filename = 'Valve_Data/ValveDirichlet.mph';
model = mphopen(filename);

A_comsol  = mphinterp(model, 'mf.Az', 'coord', mesh.nodes2coord', 'solnum', 1)';
Bx_comsol = mphinterp(model, 'mf.Bx', 'coord', mesh.nodes2coord', 'solnum', 1)';
By_comsol = mphinterp(model, 'mf.By', 'coord', mesh.nodes2coord', 'solnum', 1)';

zoom = [params.x_piston_min - tol, params.x_piston_max + tol; params.y_piston_min - tol, params.y_piston_max + tol];

%A
filename = 'Valve_Results/Dirichlet2_nonlinear/A0.png';
fieldname = 'A';
field1 = A_comsol;
field2 = A;
SubPlotData(x,y,elems2nodes,field1,field2,params,[],filename,fieldname,save);

filename = 'Valve_Results/Dirichlet2_nonlinear/A1.png';
fieldname = 'A';
field1 = A_comsol;
field2 = A;
SubPlotData(x,y,elems2nodes,field1,field2,params,zoom,filename,fieldname,save);
%Bx
filename = 'Valve_Results/Dirichlet2_nonlinear/Bx0.png';
fieldname = 'Bx';
field1 = Bx_comsol;
field2 = B(:,1);
SubPlotData(x,y,elems2nodes,field1,field2,params,[],filename,fieldname,save);

filename = 'Valve_Results/Dirichlet2_nonlinear/Bx1.png';
fieldname = 'Bx';
field1 = Bx_comsol;
field2 = B(:,1);
SubPlotData(x,y,elems2nodes,field1,field2,params,zoom,filename,fieldname,save);

%By
filename = 'Valve_Results/Dirichlet2_nonlinear/By0.png';
fieldname = 'By';
field1 = By_comsol;
field2 = B(:,2);
SubPlotData(x,y,elems2nodes,field1,field2,params,[],filename,fieldname,save);

filename = 'Valve_Results/Dirichlet2_nonlinear/By1.png';
fieldname = 'By';
field1 = By_comsol;
field2 = B(:,2);
SubPlotData(x,y,elems2nodes,field1,field2,params,zoom,filename,fieldname,save);



