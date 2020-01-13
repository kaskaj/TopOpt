function Motor_MeshGen(params, show)

if nargin < 2
    show = 0;
end

%% Modify the required edges inside

part{1} = params.edges{1}(:,1);     %Stator (upper part)
part{2} = params.edges{2}(:,1);     %Stator (lower part)
part{3} = [...                      %Rotor with hole
           params.edges{3}(:,1)
           params.edges{6}(:,1)
          ];
part{4} = params.edges{4}(:,1);     %Airgap
part{5} = params.edges{5}(:,1);     %Shaft
part{6} = params.edges{6}(:,1);     %Hole

for i = 7:length(params.edges)    %Slots
    part{i} = params.edges{i}(:,1);
end

edge = params.edges{1}(:,1:2);
node = params.nodes{1};

for i = 2:length(params.edges)
    edge = [edge;params.edges{i}(:,1:2)];    
end

for i = 2:length(params.edges)
    node = [node;params.nodes{i}];    
end
%% Add these edges and create mesh

hmax = +0.5e-1;
[vlfs,tlfs,hlfs] = lfshfn2(node,edge,part) ;

hlfs = min(hmax,hlfs);
slfs = idxtri2(vlfs,tlfs);
hfun = @trihfn2;

[nodes2coord,bedges2nodes,elems2nodes,tnum] = refine2(node,edge,part,[],hfun,vlfs,tlfs,slfs,hlfs);


%% Possibly show mesh

if show
    show_mesh(elems2nodes,nodes2coord);
    Motor_PlotEdges(params, 0, 'r', 2);
end

%% Save mesh

mesh = [];
mesh.nodes2coord  = nodes2coord;
mesh.bedges2nodes = bedges2nodes;
mesh.elems2nodes  = elems2nodes;
mesh.tnum         = tnum;

file_name = fullfile('Motor_Data', 'Mesh0.mat');
save(file_name, 'mesh');

end


