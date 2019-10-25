function [] = Valve_MeshGen(show)

addpath(genpath('./3rd_party'));

load(fullfile('Valve_Data', 'param'));

if nargin < 1
    show = 0;
end

nod1 = [
    0, -y_air/2;      % Air region
    0, y_air/2;
    x_air, y_air/2;
    x_air, -y_air/2;
    ] ;
edg1 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
edg1(:,3) = +0;

nod2 = [
    x_channel + x1 + 2*t, (y3/2) + t + y2;        % Conductor 1 (upper)
    x_channel + x1 + 2*t, (y3/2) + t;
    x_channel + x1 + 2*t + x2, (y3/2) + t;
    x_channel + x1 + 2*t + x2, (y3/2) + t + y2;
    ] ;
edg2 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
edg2(:,3) = +1;

nod3 = [
    (x_channel + x1 +2*t), -((y3/2) + t + y2);        % Conductor 2 (lower)
    (x_channel + x1 +2*t), -((y3/2) + t);
    (x_channel + x1 +2*t + x2), -((y3/2) + t);
    (x_channel + x1 +2*t + x2), -((y3/2) + t + y2);
    ] ;
edg3 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
edg3(:,3) = +2;

nod4 = [
    (x_channel),  -(y3/2 + 2*t + y2 - y4 - t);        % Plunger
    (x_channel + x1),  -(y3/2 + 2*t + y2 - y4 - t);
    (x_channel + x1),  -(y3/2 + 2*t + y2 - y4 - t - y_piston);
    (x_channel),  -(y3/2 + 2*t + y2 - y4 - t - y_piston);
    ] ;
edg4 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
edg4(:,3) = +3;


nod5 = [
    x_channel, (y3/2) + y2 + 2*t - y4;        % Iron
    x_channel, (y3/2) + y2 + 2*t + y1;
    x_channel + x1 + x4 + w_m + t + x3, (y3/2) + y2 + 2*t + y1;
    x_channel + x1 + x4 + w_m + t + x3, -((y3/2) + y2 + 2*t + y1);
    x_channel, -((y3/2) + y2 + 2*t + y1);
    x_channel, -((y3/2) + y2 + 2*t - y4);
    x_channel + x1, -((y3/2) + y2 + 2*t - y4);
    x_channel + x1, -((y3/2) + y2 + 2*t);
    x_channel + x1 + x4 + w_m + 2*t , -((y3/2) + y2 + 2*t);
    x_channel + x1 + x4 + w_m + 2*t , -((y3/2));
    x_channel + x1 + x4 + w_m - x2 , -((y3/2));
    x_channel + x1 + x4 + w_m - x2 , (y3/2);
    x_channel + x1 + x4 + w_m + 2*t , (y3/2);
    x_channel + x1 + x4 + w_m + 2*t , (y3/2) + 2*t + y2;
    x_channel + x1 , (y3/2) +  2*t + y2;
    x_channel + x1 , (y3/2) + 2*t + y2 - y4;
    ] ;
edg5 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  5
    5 ,  6 ;  6 ,  7
    7 ,  8 ;  8 ,  9
    9 ,  10 ;  10 ,  11
    11 ,  12 ;  12 ,  13
    13 ,  14 ;  14 ,  15
    15 ,  16 ;  16 ,  1
    ] ;
edg5(:,3) = +4;



edg2(:,1:2) = edg2(:,1:2) + size(nod1,1);
edg3(:,1:2) = edg3(:,1:2) + size(nod1,1) + size(nod2,1);
edg4(:,1:2) = edg4(:,1:2) + size(nod1,1) + size(nod2,1) + size(nod3,1);
edg5(:,1:2) = edg5(:,1:2) + size(nod1,1) + size(nod2,1) + size(nod3,1) + size(nod4,1);

edge = [edg1; edg2; edg3; edg4; edg5];
node = [nod1; nod2; nod3; nod4; nod5];

part{1} = [ ...
    find(edge(:,3) == 0)
    find(edge(:,3) == 1)
    find(edge(:,3) == 2)
    find(edge(:,3) == 3)
    find(edge(:,3) == 4)   
    ] ;
part{2} = [ ...
    find(edge(:,3) == 1)
    ] ;
part{3} = [ ...
    find(edge(:,3) == 2)
    ] ;
part{4} = [ ...
    find(edge(:,3) == 3)
    ];
part{5} = [ ...
    find(edge(:,3) == 4)
    ] ;


edge = edge(:,1:2) ;

hmax = +0.05;

[vlfs,tlfs,hlfs] = lfshfn2(node,edge,part) ;

hlfs = min(hmax,hlfs) ;

[slfs] = idxtri2(vlfs,tlfs) ;

%---------------------------------------------- do mesh-gen.
hfun = @trihfn2;

[nodes2coord,bedges2nodes,elems2nodes,tnum] = refine2(node,edge,part,[],hfun,vlfs,tlfs,slfs,hlfs);

if show
    show_mesh(elems2nodes,nodes2coord);    
end

file_name = fullfile('Valve_Data', 'mesh0.mat');
save(file_name,'nodes2coord','bedges2nodes','elems2nodes');
end


