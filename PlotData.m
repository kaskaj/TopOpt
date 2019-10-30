function [] = PlotData(x,y,elems2nodes,R,file_name,save)
    
    if nargin < 6
        save = 0;
    end
    
    fig1 = figure();
    trisurf(elems2nodes, x, y, R);
    view(2);
    colorbar;
    xlim([min(x), max(x)]);
    ylim([min(y), max(y)]);
    shading interp;
    colormap jet;
    axis equal;
    
    if save
        saveas(fig1, file_name);
    end
    
end

