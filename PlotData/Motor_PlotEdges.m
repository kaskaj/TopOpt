function [] = Motor_PlotEdges(params, z, colortag, lw)

if nargin < 4
    lw = 1;
end
if nargin < 3
    colortag = 'k-';
end

for i = 1:length(params.nodes)
    plot.nodes{i}  = [params.nodes{i}; params.nodes{i}(1,:)];
end

if nargin < 1
    hold on;
    
    for i = 1:length(plot.nodes)
        plot(plot.nodes{i}(:,1), plot.nodes{i}(:,2),colortag,'LineWidth',lw);        
    end
    
    axis equal
else
    hold on;
    
    for i = 1:length(plot.nodes)
        plot3(plot.nodes{i}(:,1),plot.nodes{i}(:,2),repmat(z,size(plot.nodes{i},1),1),colortag,'LineWidth',lw);       
    end
    
end
end





