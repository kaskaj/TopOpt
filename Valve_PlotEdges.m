function [] = Valve_PlotEdges(params, z, colortag, lw)

if nargin < 4
    lw = 1;
end
if nargin < 3
    colortag = 'k-';
end


nod1 = [params.nod1; params.nod1(1,:)];
nod2 = [params.nod2; params.nod2(1,:)];
nod3 = [params.nod3; params.nod3(1,:)];
nod4 = [params.nod4; params.nod4(1,:)];
nod5 = [params.nod5; params.nod5(1,:)];

if nargin < 1
    hold on;
    plot(nod1(:,1),nod1(:,2),colortag,'LineWidth',lw);
    plot(nod2(:,1),nod2(:,2),colortag,'LineWidth',lw);
    plot(nod3(:,1),nod3(:,2),colortag,'LineWidth',lw);
    plot(nod4(:,1),nod4(:,2),colortag,'LineWidth',lw);
    plot(nod5(:,1),nod5(:,2),colortag,'LineWidth',lw);
    axis equal
else
    hold on;
    plot3(nod1(:,1),nod1(:,2),repmat(z,size(nod1,1),1),colortag,'LineWidth',lw);
    plot3(nod2(:,1),nod2(:,2),repmat(z,size(nod2,1),1),colortag,'LineWidth',lw);
    plot3(nod3(:,1),nod3(:,2),repmat(z,size(nod3,1),1),colortag,'LineWidth',lw);
    plot3(nod4(:,1),nod4(:,2),repmat(z,size(nod4,1),1),colortag,'LineWidth',lw);
    plot3(nod5(:,1),nod5(:,2),repmat(z,size(nod5,1),1),colortag,'LineWidth',lw);
end
end





