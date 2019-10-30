function show_mesh(elems2nodes,nodes2coord)
if (size(elems2nodes,2)==3)
    X=reshape(nodes2coord(elems2nodes',1),size(elems2nodes,2),size(elems2nodes,1));
    Y=reshape(nodes2coord(elems2nodes',2),size(elems2nodes,2),size(elems2nodes,1));
    patch(X,Y,1);
    axis equal;
else
    tetramesh(elems2nodes,nodes2coord,'FaceAlpha',1);camorbit(20,0);
    axis equal;
end


