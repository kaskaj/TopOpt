function [J,phi] = Motor_MoveCurrent(mesh, params, step, t)

Im = params.Im;
S  = params.S;
N  = params.N;
w0 = params.w0;

Ia  = Im*cos(w0*t);                 
Ib  = Im*cos(w0*t - 2*pi/3);       
Ic  = Im*cos(w0*t + 2*pi/3);      

Ja = N*Ia*S;
Jb = N*Ib*S;
Jc = N*Ic*S;

%% Prescribe current density

J_ele = zeros(mesh.nelement,1);
J     = zeros(mesh.npoint,1);

pos = [7:43]';
pos = circshift(pos,step);

ii_a = ismember(mesh.tnum, [pos(1:3),pos(7:9)]);
ii_c = ismember(mesh.tnum, [pos(13:15),pos(19:21)]);
ii_b = ismember(mesh.tnum, [pos(25:27),pos(31:33)]);

J_ele(ii_a)  = 1;
J_ele(ii_c)  = 2;
J_ele(ii_b)  = 3;

J(mesh.elems2nodes(J_ele == 1)) =  Ja;
J(mesh.elems2nodes(J_ele == 2)) = -Jc;
J(mesh.elems2nodes(J_ele == 3)) =  Jb;


%% Prescribe fixed Air and Iron domains

phi = zeros(mesh.nelement,1);

ii_fix0 = ismember(mesh.tnum, [4,6]);
ii_fix1 = ismember(mesh.tnum, [1,2,3,5]);
ii_coil = ismember(mesh.tnum, 7:72);
ii_tooth = ismember(mesh.tnum, 7:72);

ii_coil = ii_coil & (J_ele ~= 0);
ii_tooth = ii_tooth & (J_ele == 0);

ii_fix0 = ii_fix0 | ii_coil;
ii_fix1 = ii_fix1 | ii_tooth;

ii_fix  = ii_fix0 | ii_fix1;
ii_opt  = ~ii_fix;
 
phi(ii_fix0)  = 0;
phi(ii_fix1)  = 1;

%Plot prescribed domains
% figure;
% Motor_PlotEdges(params, 1, 'k-', 2);
% plot(mesh.x_mid(ii_fix0),mesh.y_mid(ii_fix0),'o','MarkerFaceColor','b');
% plot(mesh.x_mid(ii_fix1),mesh.y_mid(ii_fix1),'ro','MarkerFaceColor','r');
% axis equal;
% 
% PlotData(mesh.x,mesh.y,mesh.elems2nodes,J);
% Motor_PlotEdges(params,max(J));

end

