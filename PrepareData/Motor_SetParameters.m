%% Set parameters

d	= 1e-3;             %Airgap
L	= 81e-3;            %Length
D1	= 72e-3;            %Stator (inner)
D2	= 110e-3;           %Stator (outer)
D3	= D1 - (2*d);       %Rotor (outer)
Dh	= 22e-3;            %Shaft
Dc  = 10e-3;            %Hole diameter

Im = 4*sqrt(2)*1e6;
w0 = 2*pi*50;
t = 0;
Ja  = Im*cos(w0*t);                   % supply current - coil 1
Jb  = Im*cos(w0*t - 2*pi/3);          % supply current - coil 2
Jc  = Im*cos(w0*t + 2*pi/3);          % supply current - coil 3
mu0 = 4*pi*1e-7;
nslots = 36;
mur = 700;
alpha = (90/nslots)*(pi/180);
gamma = pi/600;

tol1 = 3e-3;
tol2 = 1e-3;
tol3 = 0.25e-3;

params = [];
params.D1        = D1;
params.D2        = D2;
params.D3        = D3;
params.Dh        = Dh;
params.Dc        = Dc;
params.d         = d;
params.L         = L;
params.Ja        = Ja;
params.Jb        = Jb;
params.Jc        = Jc;
params.w0        = w0;
params.Im        = Im;
params.mu0       = mu0;
params.mur       = mur;
params.nslots    = nslots;

%% Compute edges

points = 50;
points0 = nslots;
points2 = points/5;

phi0 = linspace(0,pi/2,nslots+1);
phi = linspace(0,pi/2,points);
phi2 = linspace(0,pi/2,points2);
phi3 = linspace(0,pi,points2);

params.nodes  = [];
params.edges = [];

% %Stator
% points = 50;
% phi = linspace(0,pi/2,points);
% x = [(params.D2/2)*cos(phi),(params.D1/2)*cos(flip(phi))];
% y = [(params.D2/2)*sin(phi),(params.D1/2)*sin(flip(phi))];
% e1 = (1:2*points);
% e2 = [(2:2*points),1];
% params.nodes{1} = [x',y'] ;
% params.edges{1} = [e1',e2'];

%Stator (upper part)
x = [(params.D2/2)*cos(phi),(params.D2/2 - tol1)*cos(flip(phi0))];
y = [(params.D2/2)*sin(phi),(params.D2/2 - tol1)*sin(flip(phi0))];
e1 = [(1:points),(points+1):(points+1+points0)];
e2 = [(2:points),(points+1):(points+1+points0),1];
params.nodes{1} = [x',y'] ;
params.edges{1} = [e1',e2'];

%Stator (lower part)
x = [(params.D1/2 + tol2)*cos(phi0),(params.D1/2)*cos(flip(phi))];
y = [(params.D1/2 + tol2)*sin(phi0),(params.D1/2)*sin(flip(phi))];
e1 = [(1:points),(points+1):(points+1+points0)];
e2 = [(2:points),(points+1):(points+1+points0),1];
params.nodes{2} = [x',y'] ;
params.edges{2} = [e1',e2'];

%Rotor
x = [(params.D3/2)*cos(phi),(params.Dh/2)*cos(flip(phi2))];
y = [(params.D3/2)*sin(phi),(params.Dh/2)*sin(flip(phi2))];
e1 = (1:points + points2);
e2 = [(2:points + points2),1];
params.nodes{3} = [x',y'] ;
params.edges{3} = [e1',e2'];

%Airgap
x = [(params.D1/2)*cos(phi),(params.D3/2)*cos(flip(phi))];
y = [(params.D1/2)*sin(phi),(params.D3/2)*sin(flip(phi))];
e1 = (1:2*points);
e2 = [(2:2*points),1];
params.nodes{4} = [x',y'] ;
params.edges{4} = [e1',e2'];

%Shaft
x = [(params.Dh/2)*cos(phi2),zeros(1,points2)];
y = [(params.Dh/2)*sin(phi2),zeros(1,points2)];
e1 = (1:2*points2);
e2 = [(2:2*points2),1];
params.nodes{5} = [x',y'] ;
params.edges{5} = [e1',e2'];

%Hole
x0 = 15e-3;
y0 = 15e-3;
x = [x0+(params.Dc/2)*cos(phi3),x0-(params.Dc/2)*cos(phi3)];
y = [y0+(params.Dc/2)*sin(phi3),y0-(params.Dc/2)*sin(phi3)];
e1 = (1:2*points2);
e2 = [(2:2*points2),1];
params.nodes{6} = [x',y'] ;
params.edges{6} = [e1',e2'];

%Slots
R1 = params.D1/2+tol2;
R2 = params.D2/2-tol1;
e1 = (1:4);
e2 = [(2:4),1];

for i=1:length(phi0)-1
   params.nodes{6+i} = [R1*cos(phi0(i)),R1*sin(phi0(i));
                        R1*cos(phi0(i+1)),R1*sin(phi0(i+1));
                        R2*cos(phi0(i+1)),R2*sin(phi0(i+1));
                        R2*cos(phi0(i)),R2*sin(phi0(i));              
                       ]; 
   params.edges{6+i} = [e1',e2'];
end

%Slots - smaller
R1 = params.D1/2+tol2+tol3;
R2 = params.D2/2-tol1-tol3;
e1 = (1:4);
e2 = [(2:4),1];

for i=1:length(phi0)-1
   params.nodes{6+length(phi0)-1+i} = [R1*cos(phi0(i)+gamma),R1*sin(phi0(i)+gamma);
                                       R1*cos(phi0(i+1)-gamma),R1*sin(phi0(i+1)-gamma);
                                       R2*cos(phi0(i+1)-gamma),R2*sin(phi0(i+1)-gamma);
                                       R2*cos(phi0(i)+gamma),R2*sin(phi0(i)+gamma);              
                                       ]; 
   params.edges{6+length(phi0)-1+i} = [e1',e2'];
end


len = 0;
for i=1:length(params.edges)
    aux = params.edges{i};
    aux(:,1:2) = aux(:,1:2) + len;
    aux(:,3)   = i-1;
    params.edges{i} = aux;

    len = len + size(params.nodes{i},1);
end

%% Create the structure

file_name = fullfile('Motor_Data', 'Param.mat');
save(file_name, 'params');

Motor_PlotEdges(params, 0)
axis equal;

