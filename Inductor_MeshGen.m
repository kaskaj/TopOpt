function [] = Inductor_MeshGen(show)

addpath(genpath('./3rd_party'));

load(fullfile('Inductor_Data', 'param'));

if nargin < 1
    show = 0;
end

nod1 = [
    -(d1 + 2*L)/2, -d1/2;      % Air region
    -(d1 + 2*L)/2,  d1/2;
    (d1 + 2*L)/2,  d1/2;
    (d1 + 2*L)/2, -d1/2;
    ] ;
edg1 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
edg1(:,3) = +0;

nod2 = [
    -2*L - d2/2, -d2/2;        % Conductor 1 (left)
    -2*L - d2/2,  d2/2;
    -2*L + d2/2,  d2/2;
    -2*L + d2/2, -d2/2;
    ] ;
edg2 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
edg2(:,3) = +1;

nod3 = [
    2*L - d2/2, -d2/2;        % Conductor 2 (right)
    2*L - d2/2,  d2/2;
    2*L + d2/2,  d2/2;
    2*L + d2/2, -d2/2;
    ] ;
edg3 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
edg3(:,3) = +2;

nod4 = [
    -d3/2,  -h3/2;        % Iron
    -d3/2,   h3/2;
    d3/2,   h3/2;
    d3/2,  -h3/2;
    ] ;
edg4 = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ] ;
edg4(:,3) = +3;

edg2(:,1:2) = edg2(:,1:2) + size(nod1,1);
edg3(:,1:2) = edg3(:,1:2) + size(nod1,1) + size(nod2,1);
edg4(:,1:2) = edg4(:,1:2) + size(nod1,1) + size(nod2,1) + + size(nod3,1);

edge = [edg1; edg2; edg3; edg4];
node = [nod1; nod2; nod3; nod4];

part{1} = [ ...
    find(edge(:,3) == 0)
    find(edge(:,3) == 1)
    find(edge(:,3) == 2)
    find(edge(:,3) == 3)
    ] ;
part{2} = [ ...
    find(edge(:,3) == 1)
    ] ;
part{3} = [ ...
    find(edge(:,3) == 2)
    ] ;
part{4} = [ ...
    find(edge(:,3) == 3)
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

file_name = fullfile('Inductor_Data', 'mesh0.mat');
save(file_name,'nodes2coord','bedges2nodes','elems2nodes');
end


