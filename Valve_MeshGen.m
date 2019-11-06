function Valve_MeshGen(params, show)
    
    addpath(genpath('./3rd_party'));
    
    if nargin < 2
        show = 0;
    end
    
    %% Modify the required edges inside
    
    edge = [params.edg1; params.edg2; params.edg3; params.edg4; params.edg5; params.edg6; params.edg7; params.edg8;];
    node = [params.nod1; params.nod2; params.nod3; params.nod4; params.nod5; params.nod6; params.nod7; params.nod8;];
    
    part{1} = [ ...
        find(edge(:,3) == 0)
        find(edge(:,3) == 1)
        find(edge(:,3) == 2)
        find(edge(:,3) == 3)
        find(edge(:,3) == 4)
        find(edge(:,3) == 5)
        find(edge(:,3) == 6)
        find(edge(:,3) == 7)
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
     part{6} = [ ...
        find(edge(:,3) == 5)
        ] ;
     part{7} = [ ...
        find(edge(:,3) == 6)
        ] ;
     part{8} = [ ...
        find(edge(:,3) == 7)
        ] ;
    
    edge = edge(:,1:2) ;
    
    %% Add these edges and create mesh
    
    hmax = +0.05;
    
    [vlfs,tlfs,hlfs] = lfshfn2(node,edge,part) ;
    
    hlfs = min(hmax,hlfs);
    slfs = idxtri2(vlfs,tlfs);
    hfun = @trihfn2;
    
    [nodes2coord,bedges2nodes,elems2nodes,tnum] = refine2(node,edge,part,[],hfun,vlfs,tlfs,slfs,hlfs);
    
    %% Possibly show mesh
    
    if show
        show_mesh(elems2nodes,nodes2coord);
    end
    
    %% Save mesh
    
    mesh = [];
    mesh.nodes2coord  = nodes2coord;
    mesh.bedges2nodes = bedges2nodes;
    mesh.elems2nodes  = elems2nodes;
    mesh.tnum         = tnum;
    
    file_name = fullfile('Valve_Data', 'Mesh0.mat');
    save(file_name, 'mesh');
    
end


