function [] = SubPlotMesh(elems2nodes1,nodes2coord1,elems2nodes2,nodes2coord2,params,zoom,file_name,mesh_name,save)

if isempty(zoom)
    zoom = [min(nodes2coord1(:,1)), max(nodes2coord1(:,1)); min(nodes2coord1(:,2)), max(nodes2coord1(:,2))];
end
if nargin < 9
    save = 0;
end

fig1 = figure();
subplot(1,2,1)
X=reshape(nodes2coord1(elems2nodes1',1),size(elems2nodes1,2),size(elems2nodes1,1));
Y=reshape(nodes2coord1(elems2nodes1',2),size(elems2nodes1,2),size(elems2nodes1,1));
patch(X,Y,1);
axis equal;
Valve_PlotEdges(params, 0, 'r-', 2);
title(strcat('Comsol: ',mesh_name))
xlim(zoom(1,:));
ylim(zoom(2,:));

subplot(1,2,2)
X=reshape(nodes2coord2(elems2nodes2',1),size(elems2nodes2,2),size(elems2nodes2,1));
Y=reshape(nodes2coord2(elems2nodes2',2),size(elems2nodes2,2),size(elems2nodes2,1));
patch(X,Y,1);
axis equal;
Valve_PlotEdges(params, 0, 'r-', 2);
title(strcat('Matlab: ',mesh_name))
xlim(zoom(1,:));
ylim(zoom(2,:));

if save
    saveas(fig1, file_name);
end

end