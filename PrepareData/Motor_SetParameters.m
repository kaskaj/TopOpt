%% Set parameters

d	= 1e-3;             %Airgap
L	= 81e-3;             %Length
D1	= 72e-3;            %Stator (inner)
D2	= 110e-3;           %Stator (outer)
D3	= D1 - (2*d);    %Rotor (outer)
Dh	= 22e-3;            %Shaft

Im = 4*sqrt(2);
w0 = 2*pi*50;
t = 0;
Ia  = Im*cos(w0*t);                   % supply current - coil 1
Ib  = Im*cos(w0*t - 2*pi/3);          % supply current - coil 2
Ic  = Im*cos(w0*t + 2*pi/3);          % supply current - coil 3
N   = 470;                            % number of turns
mu0 = 4*pi*1e-7;
mur = 700;
S = 1.0271e2;                        %Slot area

Ja = (N*Ia)*S;
Jb = (N*Ib)*S;
Jc = (N*Ic)*S;

params = [];
params.D1        = D1;
params.D2        = D2;
params.D3        = D3;
params.Dh        = Dh;
params.d         = d;
params.L         = L;
params.Ja        = Ja;
params.Jb        = Jb;
params.Jc        = Jc;
params.mu0       = mu0;
params.mur       = mur;

%% Compute edges

%Stator
points = 50;
phi = linspace(0,pi/2,points);
x = [(params.D2/2)*cos(phi),(params.D1/2)*cos(flip(phi))];
y = [(params.D2/2)*sin(phi),(params.D1/2)*sin(flip(phi))];
params.nod1 = [x',y'] ;
e1 = (1:2*points);
e2 = [(2:2*points),1];
params.edg1 = [e1',e2'];
params.edg1(:,3) = +0;

%Slot1
alpha = pi/18;
beta = pi/18/4;
tol = 3e-3;
phi0 = -(pi/2)/6 + beta;
params.nod2 = [(params.D1/2+tol)*cos(phi0 + pi/2),(params.D1/2+tol)*sin(phi0 + pi/2);
               (params.D1/2+tol)*cos(phi0 + pi/2 + alpha),(params.D1/2+tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2 + alpha),(params.D2/2-tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2),(params.D2/2-tol)*sin(phi0 + pi/2);              
               ] ;
e1 = (1:4);
e2 = [(2:4),1];
params.edg2 = [e1',e2'];
params.edg2(:,3) = +1;

%Slot2
phi0 = -2*(pi/2)/6 + beta;
params.nod3 = [(params.D1/2+tol)*cos(phi0 + pi/2),(params.D1/2+tol)*sin(phi0 + pi/2);
               (params.D1/2+tol)*cos(phi0 + pi/2 + alpha),(params.D1/2+tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2 + alpha),(params.D2/2-tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2),(params.D2/2-tol)*sin(phi0 + pi/2);              
               ] ;
e1 = (1:4);
e2 = [(2:4),1];
params.edg3 = [e1',e2'];
params.edg3(:,3) = +2;

%Slot3
phi0 = -3*(pi/2)/6 + beta;
params.nod4 = [(params.D1/2+tol)*cos(phi0 + pi/2),(params.D1/2+tol)*sin(phi0 + pi/2);
               (params.D1/2+tol)*cos(phi0 + pi/2 + alpha),(params.D1/2+tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2 + alpha),(params.D2/2-tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2),(params.D2/2-tol)*sin(phi0 + pi/2);              
               ] ;
e1 = (1:4);
e2 = [(2:4),1];
params.edg4 = [e1',e2'];
params.edg4(:,3) = +3;

%Slot4
phi0 = -4*(pi/2)/6 + beta;
params.nod5 = [(params.D1/2+tol)*cos(phi0 + pi/2),(params.D1/2+tol)*sin(phi0 + pi/2);
               (params.D1/2+tol)*cos(phi0 + pi/2 + alpha),(params.D1/2+tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2 + alpha),(params.D2/2-tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2),(params.D2/2-tol)*sin(phi0 + pi/2);              
               ] ;
e1 = (1:4);
e2 = [(2:4),1];
params.edg5 = [e1',e2'];
params.edg5(:,3) = +4;

%Slot5
phi0 = -5*(pi/2)/6 + beta;
params.nod6 = [(params.D1/2+tol)*cos(phi0 + pi/2),(params.D1/2+tol)*sin(phi0 + pi/2);
               (params.D1/2+tol)*cos(phi0 + pi/2 + alpha),(params.D1/2+tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2 + alpha),(params.D2/2-tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2),(params.D2/2-tol)*sin(phi0 + pi/2);              
               ] ;
e1 = (1:4);
e2 = [(2:4),1];
params.edg6 = [e1',e2'];
params.edg6(:,3) = +5;

%Slot6
phi0 = -6*(pi/2)/6 + beta;
params.nod7 = [(params.D1/2+tol)*cos(phi0 + pi/2),(params.D1/2+tol)*sin(phi0 + pi/2);
               (params.D1/2+tol)*cos(phi0 + pi/2 + alpha),(params.D1/2+tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2 + alpha),(params.D2/2-tol)*sin(phi0 + pi/2 + alpha);
               (params.D2/2-tol)*cos(phi0 + pi/2),(params.D2/2-tol)*sin(phi0 + pi/2);              
               ] ;
e1 = (1:4);
e2 = [(2:4),1];
params.edg7 = [e1',e2'];
params.edg7(:,3) = +6;

%Airgap
points = 50;
phi = linspace(0,pi/2,points);
x = [(params.D1/2)*cos(phi),(params.D3/2)*cos(flip(phi))];
y = [(params.D1/2)*sin(phi),(params.D3/2)*sin(flip(phi))];
params.nod8 = [x',y'] ;
e1 = (1:2*points);
e2 = [(2:2*points),1];
params.edg8 = [e1',e2'];
params.edg8(:,3) = +7;

%Rotor
points = 50;
phi = linspace(0,pi/2,points);
phi2 = linspace(0,pi/2,points/5);
x = [(params.D3/2)*cos(phi),(params.Dh/2)*cos(flip(phi2))];
y = [(params.D3/2)*sin(phi),(params.Dh/2)*sin(flip(phi2))];
params.nod9 = [x',y'] ;
e1 = (1:points + points/5);
e2 = [(2:points + points/5),1];
params.edg9 = [e1',e2'];
params.edg9(:,3) = +8;

%Shaft
points = 10;
phi = linspace(0,pi/2,points);
x = [(params.Dh/2)*cos(phi),zeros(1,length(phi))];
y = [(params.Dh/2)*sin(phi),zeros(1,length(phi))];
params.nod10 = [x',y'] ;
e1 = (1:2*points);
e2 = [(2:2*points),1];
params.edg10 = [e1',e2'];
params.edg10(:,3) = +9;

params.edg2(:,1:2 ) = params.edg2(:,1:2) + size(params.nod1,1);
params.edg3(:,1:2)  = params.edg3(:,1:2) + size(params.nod1,1) + size(params.nod2,1);
params.edg4(:,1:2)  = params.edg4(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1);
params.edg5(:,1:2)  = params.edg5(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1);
params.edg6(:,1:2)  = params.edg6(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1) + size(params.nod5,1);
params.edg7(:,1:2)  = params.edg7(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1) + size(params.nod5,1) + size(params.nod6,1);
params.edg8(:,1:2)  = params.edg8(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1) + size(params.nod5,1) + size(params.nod6,1) + size(params.nod7,1);
params.edg9(:,1:2)  = params.edg9(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1) + size(params.nod5,1) + size(params.nod6,1) + size(params.nod7,1) + size(params.nod8,1);
params.edg10(:,1:2) = params.edg10(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1) + size(params.nod5,1) + size(params.nod6,1) + size(params.nod7,1) + size(params.nod8,1) + size(params.nod9,1);
 

%% Create the structure

file_name = fullfile('Motor_Data', 'Param.mat');
save(file_name, 'params');


