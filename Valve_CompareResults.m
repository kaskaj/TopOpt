load(fullfile('Valve_Data', 'Mesh5'));
load(fullfile('Valve_Data', 'Param'))

filename = 'Valve_Data/Valve.mph';
model = mphopen(filename);

A_comsol  = mphinterp(model, 'mf.Az', 'coord', mesh.nodes2coord', 'solnum', 1)';
Bx_comsol = mphinterp(model, 'mf.Bx', 'coord', mesh.nodes2coord', 'solnum', 1)';
By_comsol = mphinterp(model, 'mf.By', 'coord', mesh.nodes2coord', 'solnum', 1)';


tol = 0.001;
zoom = [params.x_piston_min - tol, params.x_piston_max + tol; params.y_piston_min - tol, params.y_piston_max + tol];

%A
filename = 'Valve_Results/A0.png';
fieldname = 'A';
field1 = A_comsol;
field2 = A;
SubPlotData(x,y,elems2nodes,field1,field2,params,[],filename,fieldname,1);

filename = 'Valve_Results/A1.png';
fieldname = 'A';
field1 = A_comsol;
field2 = A;
SubPlotData(x,y,elems2nodes,field1,field2,params,zoom,filename,fieldname,1);
%Bx
filename = 'Valve_Results/Bx0.png';
fieldname = 'Bx';
field1 = Bx_comsol;
field2 = B(:,1);
SubPlotData(x,y,elems2nodes,field1,field2,params,[],filename,fieldname,1);

filename = 'Valve_Results/Bx1.png';
fieldname = 'Bx';
field1 = Bx_comsol;
field2 = B(:,1);
SubPlotData(x,y,elems2nodes,field1,field2,params,zoom,filename,fieldname,1);

%By
filename = 'Valve_Results/By0.png';
fieldname = 'By';
field1 = By_comsol;
field2 = B(:,2);
SubPlotData(x,y,elems2nodes,field1,field2,params,[],filename,fieldname,1);

filename = 'Valve_Results/By1.png';
fieldname = 'By';
field1 = By_comsol;
field2 = B(:,2);
SubPlotData(x,y,elems2nodes,field1,field2,params,zoom,filename,fieldname,1);



