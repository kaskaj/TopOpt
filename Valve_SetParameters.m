%% Set parameters

I1  = 2;         % supply current - coil 1
I2  = 2;         % supply current - coil 2
N   = 900;       % number of turns
t   = 0.001;     % technological size
mu0 = 4*pi*1e-7;
mur = 700;

w_magnet  = 0.0025;   % width of magnet
w_channel = 0.01;     % diameter of fluid pipe
w1 = 0.0086;          % width of inner MO part
w2 = 0.0135;          % width of coil
w3 = 0.005;           % width of outer MO part
w4 = 0.012;           % placement of magnet
h1 = 0.02;            % height of head
h2 = 0.0246;          % height of coil
h3 = 0.0104;          % height of middle MO part
h4 = 0.0210;          % height of inner MO part

h_gap    = 0.008;                          % height of air gap
h_piston = 2*h2 + h3 + 2*t - h_gap - 2*h4; % height of piston
% y_move = h_piston/4;                       % displacement of piston
y_move   = 0;

params = [];
params.I1        = I1;
params.I2        = I2;
params.N         = N;
params.t         = t;
params.mu0       = mu0;
params.mur       = mur;
params.w_magnet  = w_magnet;
params.w_channel = w_channel;
params.w1        = w1;
params.w2        = w2;
params.w3        = w3;
params.w4        = w4;
params.h1        = h1;
params.h2        = h2;
params.h3        = h3;
params.h4        = h4;
params.h_gap     = h_gap;
params.h_piston  = h_piston;


%% Compute the volumes and the current density

params.Vc = h2*pi*w2*(2*(w_channel + w1 + 2*t) + w2);            % volume of coil
params.Vm = h1*pi*w_magnet*(2*(w_channel + w1 + w4) + w_magnet); % volume of magnet

params.J_coil1 = (N*I1)/(w2*h2);
params.J_coil2 = (N*I2)/(w2*h2);


%% Compute positions of various objects

params.x_air = 8*w4;
params.y_air = 20*w4;

params.x_fe_min = w_channel;
params.x_fe_max = w_channel + w1 + 3*t + w2 + w3;
params.y_fe_min = -(h3/2 +2*t + h2 + h1);
params.y_fe_max = h3/2 +2*t + h2 + h1;

params.x_c_min = w_channel + w1;
params.x_c_max = w_channel + w1 + 3*t + w2;
params.y_c_min = h3/2;
params.y_c_max = h3/2 + 2*t + h2;

params.x_gap_min = w_channel;
params.x_gap_max = w_channel + w1;
params.y_gap_min = -(h3/2 + 2*t + h2 - h4 - t) + h_piston + y_move;
params.y_gap_max = h3/2 + 2*t + h2 - h4 - t;

params.x_t1_min = w_channel;
params.x_t1_max = w_channel + w1;
params.y_t1_min = -(h3/2 + 2*t + h2 - h4);
params.y_t1_max = -(h3/2 + 2*t + h2 - h4) + t + y_move;

params.x_t2_min = w_channel + w1;
params.x_t2_max = w_channel + w1 + t;
params.y_t2_min = -h3/2;
params.y_t2_max = h3/2;

params.x_piston_min = w_channel;
params.x_piston_max = w_channel + w1;
params.y_piston_min = -h3/2-2*t-h2+h4+t+y_move;
params.y_piston_max = -h3/2-2*t-h2+h4+t+h_piston+y_move;

%Initial iron position
params.x_init_min = w_channel;
params.x_init_max = w_channel + w1;
params.y_init_min = h3/2 + 2*t + h2 - h4;
params.y_init_max = h3/2 + 2*t + h2;


%% Compute edges

% Air region
params.nod1 = [
    0, -params.y_air/2;
    0, params.y_air/2;
    params.x_air, params.y_air/2;
    params.x_air, -params.y_air/2;
    ] ;
params.edg1 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
params.edg1(:,3) = +0;

% Coil upper
params.nod2 = [
    w_channel + w1 + 2*t, (h3/2) + t + h2;
    w_channel + w1 + 2*t, (h3/2) + t;
    w_channel + w1 + 2*t + w2, (h3/2) + t;
    w_channel + w1 + 2*t + w2, (h3/2) + t + h2;
    ] ;
params.edg2 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
params.edg2(:,3) = +1;

% Coil lower
params.nod3 = [
    (w_channel + w1 +2*t), -((h3/2) + t + h2);
    (w_channel + w1 +2*t), -((h3/2) + t);
    (w_channel + w1 +2*t + w2), -((h3/2) + t);
    (w_channel + w1 +2*t + w2), -((h3/2) + t + h2);
    ] ;
params.edg3 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
params.edg3(:,3) = +2;

% Plunger
params.nod4 = [
    (w_channel),  -(h3/2 + 2*t + h2 - h4 - t) + y_move;
    (w_channel + w1),  -(h3/2 + 2*t + h2 - h4 - t) + y_move;
    (w_channel + w1),  -(h3/2 + 2*t + h2 - h4 - t - h_piston) + y_move;
    (w_channel),  -(h3/2 + 2*t + h2 - h4 - t - h_piston) + y_move;
    ] ;
params.edg4 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
params.edg4(:,3) = +3;

% Iron
% params.nod5 = [
%     w_channel, (h3/2) + h2 + 2*t - h4;
%     w_channel, (h3/2) + h2 + 2*t + h1;
%     w_channel + w1 + 2*t + w2 + t + w3, (h3/2) + h2 + 2*t + h1;
%     w_channel + w1 + 2*t + w2 + t + w3, -((h3/2) + h2 + 2*t + h1);
%     w_channel, -((h3/2) + h2 + 2*t + h1);
%     w_channel, -((h3/2) + h2 + 2*t - h4);
%     w_channel + w1, -((h3/2) + h2 + 2*t - h4);
%     w_channel + w1, -((h3/2) + h2 + 2*t);
%     w_channel + w1 + w4 + w_magnet + 2*t , -((h3/2) + h2 + 2*t);
%     w_channel + w1 + w4 + w_magnet + 2*t , -((h3/2));
%     w_channel + w1 + w4 + w_magnet - w2 , -((h3/2));
%     w_channel + w1 + w4 + w_magnet - w2 , (h3/2);
%     w_channel + w1 + w4 + w_magnet + 2*t , (h3/2);
%     w_channel + w1 + w4 + w_magnet + 2*t , (h3/2) + 2*t + h2;
%     w_channel + w1 , (h3/2) +  2*t + h2;
%     w_channel + w1 , (h3/2) + 2*t + h2 - h4;
%     ] ;
% params.edg5 = [
%     1 ,  2 ;  2 ,  3
%     3 ,  4 ;  4 ,  5
%     5 ,  6 ;  6 ,  7
%     7 ,  8 ;  8 ,  9
%     9 ,  10 ;  10 ,  11
%     11 ,  12 ;  12 ,  13
%     13 ,  14 ;  14 ,  15
%     15 ,  16 ;  16 ,  1
%     ] ;
% params.edg5(:,3) = +4;

%Technological gap
params.nod5 = [
    w_channel + w1, h3/2 + 2*t + h2;
    w_channel + w1 + 3*t + w2, h3/2 + 2*t + h2;
    w_channel + w1 + 3*t + w2, h3/2;
    w_channel + w1 + t, h3/2;
    w_channel + w1 + t, -h3/2;
    w_channel + w1 + 3*t + w2, -h3/2;
    w_channel + w1 + 3*t + w2, -(h3/2 + 2*t + h2);
    w_channel + w1, -(h3/2 + 2*t + h2);
    %- Coil upper
    w_channel + w1 + 2*t, (h3/2) + t + h2;
    w_channel + w1 + 2*t, (h3/2) + t;
    w_channel + w1 + 2*t + w2, (h3/2) + t;
    w_channel + w1 + 2*t + w2, (h3/2) + t + h2;
    %- Coil lower
    (w_channel + w1 +2*t), -((h3/2) + t + h2);
    (w_channel + w1 +2*t), -((h3/2) + t);
    (w_channel + w1 +2*t + w2), -((h3/2) + t);
    (w_channel + w1 +2*t + w2), -((h3/2) + t + h2);
    ] ;
params.edg5 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  5
    5 ,  6 ;  6 ,  7
    7 ,  8 ;  8 ,  1
    %- Coil upper
    9 ,  10 ;  10 ,  11
    11 ,  12 ;  12 ,  9
    %- Coil lower
    13 ,  14 ;  14 ,  15
    15 ,  16 ;  16 ,  13
    ] ;
params.edg5(:,3) = +4;

% Region below the plunger
params.nod6 = [
    (w_channel),  -(h3/2 + 2*t + h2 - h4);
    (w_channel + w1),  -(h3/2 + 2*t + h2 - h4);
    (w_channel + w1),  -(h3/2 + 2*t + h2 - h4 - t - y_move);
    (w_channel),  -(h3/2 + 2*t + h2 - h4 - t - y_move);
    ] ;
params.edg6 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
params.edg6(:,3) = +5;

% Region above the plunger
params.nod7 = [
    (w_channel),  -(h3/2 + 2*t + h2 - h4 - t - y_move - h_piston);
    (w_channel + w1),  -(h3/2 + 2*t + h2 - h4 - t - y_move - h_piston);
    (w_channel + w1),  -(h3/2 + 2*t + h2 - h4 - t  - h_piston - t - h_gap);
    (w_channel),  -(h3/2 + 2*t + h2 - h4 - t  - h_piston - t - h_gap);
    ] ;
params.edg7 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
params.edg7(:,3) = +6;

% Channel
params.nod8 = [
    0,  -params.y_air/2;
    w_channel,  -params.y_air/2;
    w_channel,  params.y_air/2;
    0,  params.y_air/2;
    ] ;
params.edg8 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
params.edg8(:,3) = +7;

params.edg2(:,1:2) = params.edg2(:,1:2) + size(params.nod1,1);
params.edg3(:,1:2) = params.edg3(:,1:2) + size(params.nod1,1) + size(params.nod2,1);
params.edg4(:,1:2) = params.edg4(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1);
params.edg5(:,1:2) = params.edg5(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1);
params.edg6(:,1:2) = params.edg6(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1) + size(params.nod5,1);
params.edg7(:,1:2) = params.edg7(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1) + size(params.nod5,1) + size(params.nod6,1);
params.edg8(:,1:2) = params.edg8(:,1:2) + size(params.nod1,1) + size(params.nod2,1) + size(params.nod3,1) + size(params.nod4,1) + size(params.nod5,1) + size(params.nod6,1) + size(params.nod7,1);

%% Create the structure

file_name = fullfile('Valve_Data', 'Param.mat');
save(file_name, 'params');


