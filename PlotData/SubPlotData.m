function [] = SubPlotData(x,y,elems2nodes,R1,R2,params,zoom,file_name,field_name,save)

if isempty(zoom)
    zoom = [min(x), max(x); min(y), max(y)];
end
if nargin < 10
    save = 0;
end

ii = x >= zoom(1,1) & x <= zoom(1,2) & y >= zoom(2,1) & y <= zoom(2,2);
min_R = min([R1(ii); R2(ii)]);
max_R = max([R1(ii); R2(ii)]);

fig1 = figure();
subplot(1,2,1)
trisurf(elems2nodes, x, y, R1);
view(2);
colorbar;
shading interp;
colormap jet;
axis equal;
caxis([min_R; max_R]);
Valve_PlotEdges(params, max_R);
title(strcat('Comsol: ',field_name))
xlim(zoom(1,:));
ylim(zoom(2,:));

subplot(1,2,2)
trisurf(elems2nodes, x, y, R2)
view(2);
colorbar;
shading interp;
colormap jet;
axis equal;
caxis([min_R; max_R]);
Valve_PlotEdges(params, max_R);
title(strcat('Matlab: ',field_name))
xlim(zoom(1,:));
ylim(zoom(2,:));

if save
    saveas(fig1, file_name);
end

end