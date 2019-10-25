
cx = 0.0075; %chamfer
cy = 0.0175; %chamfer
I1 = 2; %supply current - coil 1
I2 = 2; %supply current - coil 2
N = 900; %number of turns
t = 0.001; %technological size
w_m = 0.0025; %width of magnet
x_channel = 0.01; %diameter of fluid pipe
x1 = 0.0086; %width of inner MO part
x2 = 0.0135; %width of coil
x3 = 0.005; %width of outer MO part
x4 = 0.012; %placement of magnet
y1 = 0.02; %heigth of head
y2 = 0.0246; %heigth of coil
y3 = 0.0104; %heigth of middle MO part
y4 = 0.0210; %heigth of inner MO part
y_gap = 0.008; %length of air gap
y_gap2 =  0.00795;
y_gap1 = y_gap - y_gap2;
y_piston = 2*y2 + y3 + 2*t - y_gap - 2*y4; %height of piston

Vc = y2*pi*x2*(2*(x_channel + x1 + 2*t) + x2); %Volume of coil
Vm = y1*pi*w_m*(2*(x_channel + x1 + x4) + w_m); %Volume of magnet

x_air = 8*x4;
y_air = 20*x4;

x_fe_min = x_channel;
x_fe_max = x_channel + x1 + 3*t + x2 + x3;
y_fe_min = -(y3/2 +2*t + y2 + y1);
y_fe_max = y3/2 +2*t + y2 + y1;

x_c_min = x_channel + x1;
x_c_max = x_channel + x1 + 3*t + x2;
y_c_min = y3/2;
y_c_max = y3/2 + 2*t + y2;

x_gap_min = x_channel;
x_gap_max = x_channel + x1;
y_gap_min = -(y3/2 + 2*t + y2 - y4 - t) + y_piston;
y_gap_max = y3/2 + 2*t + y2 - y4 - t;

x_t_min = x_channel;
x_t_max = x_channel + x1;
y_t_min = -(y3/2 + 2*t + y2 - y4);
y_t_max = -(y3/2 + 2*t + y2 - y4) + t;

x_t2_min = x_channel + x1;
x_t2_max = x_channel + x1 + t;
y_t2_min = -y3/2;
y_t2_max = y3/2;

mu0 = 4*pi*1e-7;
mur = 700;

J_coil1 = (N*I1)/(x2*y2);
J_coil2 = (N*I2)/(x2*y2);

file_name = fullfile('Valve_Data', 'param.mat');
save(file_name,'x1','y1','x2','y2','x3','y3','x4','y4','y_gap',...
    'y_piston','x_channel','w_m','x_air','y_air','t','mur','mu0',...
    'x_fe_min','x_fe_max','y_fe_min','y_fe_max',...
    'x_c_min','x_c_max','y_c_min','y_c_max',...
    'x_gap_min','x_gap_max','y_gap_min','y_gap_max',...
    'x_t_min','x_t_max','y_t_min','y_t_max',...
    'x_t2_min','x_t2_max','y_t2_min','y_t2_max',...
    'J_coil1','J_coil2');