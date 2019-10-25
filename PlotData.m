function [] = PlotData(x,y,R,file_name,save)

if nargin < 5
    save = 0;
end

gridDelaunay = delaunay(x,y);
fig1 = figure();
trimesh(gridDelaunay,x,y,R);
view(2);
colorbar;
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);
shading interp;
colormap jet;

if save
    saveas(fig1, file_name);
end

end

