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

coil = zeros(2*36,1);
coil(1:3)   =  1;   coil(7:9)   =  1;
coil(13:15) = -2;	coil(19:21) = -2;
coil(25:27) =  3;   coil(31:33) =  3;
coil(37:39) = -1;   coil(43:45) = -1;
coil(49:51) =  2;	coil(55:57) =  2;
coil(61:63) = -3;   coil(67:69) = -3;

pos   = (7:42);
coil   = circshift(coil,step);

for i = 1:36
   ii = ismember(mesh.tnum, pos(i));
   J_ele(ii) = coil(i);  
end

J(mesh.elems2nodes(J_ele ==  1)) =  Ja;
J(mesh.elems2nodes(J_ele == -1)) = -Ja;
J(mesh.elems2nodes(J_ele ==  2)) =  Jc;
J(mesh.elems2nodes(J_ele == -2)) = -Jc;
J(mesh.elems2nodes(J_ele ==  3)) =  Jb;
J(mesh.elems2nodes(J_ele == -3)) = -Jb;

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

coil = zeros(2*36,1);
coil(1:3)   =  1;   coil(7:9)   =  1;
coil(13:15) = -2;	coil(19:21) = -2;
coil(25:27) =  3;   coil(31:33) =  3;
coil(37:39) = -1;   coil(43:45) = -1;
coil(49:51) =  2;	coil(55:57) =  2;
coil(61:63) = -3;   coil(67:69) = -3;

pos   = (7:42);
coil   = circshift(coil,step);

for i = 1:36
   ii = ismember(mesh.tnum, pos(i));
   J_ele(ii) = coil(i);
end

J(mesh.elems2nodes(J_ele ==  1)) =  Ja;
J(mesh.elems2nodes(J_ele == -1)) = -Ja;
J(mesh.elems2nodes(J_ele ==  2)) =  Jc;
J(mesh.elems2nodes(J_ele == -2)) = -Jc;
J(mesh.elems2nodes(J_ele ==  3)) =  Jb;
J(mesh.elems2nodes(J_ele == -3)) = -Jb;

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

end

