
d1 = 1000e-3;              %Air hight
d2 = 50e-3;                %Length of the side of the condutor (square)
d3 = 150e-3;        %Width of iron
h3 = 500e-3;        %Height of iron
L = 100e-3;                %Distance of  the center of the conductor from [0,0]
I = 0.5e3;              %Current
S = d2^2;               %Surface of conductor
J_coil = I/S;       %Current density in the conductor
mu0 = 4*pi*1e-7;     %
mur = 0.5e6;
gamma_fe = 9e6;      %Electric conductivity
gamma0 = 10e-9;
f0 = 50;            %Freqency

file_name = fullfile('Inductor_Data', 'param.mat');
save(file_name,'d1','d2','d3','h3','L','I','S','J_coil','mu0',...
    'mur','gamma_fe','gamma0','f0');