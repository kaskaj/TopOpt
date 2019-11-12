function [] = Valve_PlotEdges(params, z, colortag, lw)

if nargin < 4
    lw = 1;
end
if nargin < 3
    colortag = 'k-';
end


nod1  = [params.nod1; params.nod1(1,:)];
nod2  = [params.nod2; params.nod2(1,:)];
nod3  = [params.nod3; params.nod3(1,:)];
nod4  = [params.nod4; params.nod4(1,:)];
nod5  = params.nod5;
nod6  = [params.nod6; params.nod6(1,:)];
nod7  = [params.nod7; params.nod7(1,:)];
nod8  = [params.nod8; params.nod8(1,:)];
nod9  = [params.nod9; params.nod9(1,:)];
nod10 = [params.nod10; params.nod10(1,:)];
nod11 = [params.nod11; params.nod11(1,:)];
nod12 = [params.nod12; params.nod12(1,:)];
nod13 = [params.nod13; params.nod13(1,:)];
nod14 = [params.nod14; params.nod14(1,:)];
nod15 = [params.nod15; params.nod15(1,:)];

edg5 = params.edg5;
edg9 = params.edg9;
edg10 = params.edg10;

if nargin < 1
    hold on;
    plot(nod1(:,1),nod1(:,2),colortag,'LineWidth',lw);
    plot(nod2(:,1),nod2(:,2),colortag,'LineWidth',lw);
    plot(nod3(:,1),nod3(:,2),colortag,'LineWidth',lw);
    plot(nod4(:,1),nod4(:,2),colortag,'LineWidth',lw);
    plot(nod5(:,1),nod5(:,2),colortag,'LineWidth',lw);
    plot(nod6(:,1),nod6(:,2),colortag,'LineWidth',lw);
    plot(nod7(:,1),nod7(:,2),colortag,'LineWidth',lw);
    plot(nod8(:,1),nod8(:,2),colortag,'LineWidth',lw);
    plot(nod9(:,1),nod9(:,2),colortag,'LineWidth',lw);
    plot(nod10(:,1),nod10(:,2),colortag,'LineWidth',lw);
    plot(nod11(:,1),nod11(:,2),colortag,'LineWidth',lw);
    plot(nod12(:,1),nod12(:,2),colortag,'LineWidth',lw);
    plot(nod13(:,1),nod13(:,2),colortag,'LineWidth',lw);
    plot(nod14(:,1),nod14(:,2),colortag,'LineWidth',lw);
    plot(nod15(:,1),nod15(:,2),colortag,'LineWidth',lw);
    axis equal
else
    hold on;
    plot3(nod1(:,1),nod1(:,2),repmat(z,size(nod1,1),1),colortag,'LineWidth',lw);
    plot3(nod2(:,1),nod2(:,2),repmat(z,size(nod2,1),1),colortag,'LineWidth',lw);
    plot3(nod3(:,1),nod3(:,2),repmat(z,size(nod3,1),1),colortag,'LineWidth',lw);
    plot3(nod4(:,1),nod4(:,2),repmat(z,size(nod4,1),1),colortag,'LineWidth',lw);
    for i=1:size(edg5)
        ii = edg5(i,[1 2]) - min(reshape(edg5(:,[1 2]), [], 1)) + 1;
        plot3(nod5(ii,1),nod5(ii,2),[z;z],colortag,'LineWidth',lw);
    end
    plot3(nod6(:,1),nod6(:,2),repmat(z,size(nod6,1),1),colortag,'LineWidth',lw);
    plot3(nod7(:,1),nod7(:,2),repmat(z,size(nod7,1),1),colortag,'LineWidth',lw);
    plot3(nod8(:,1),nod8(:,2),repmat(z,size(nod8,1),1),colortag,'LineWidth',lw);
    for i=1:size(edg9)
        ii = edg9(i,[1 2]) - min(reshape(edg9(:,[1 2]), [], 1)) + 1;
        plot3(nod9(ii,1),nod9(ii,2),[z;z],colortag,'LineWidth',lw);
    end
    for i=1:size(edg10)
        ii = edg10(i,[1 2]) - min(reshape(edg10(:,[1 2]), [], 1)) + 1;
        plot3(nod10(ii,1),nod10(ii,2),[z;z],colortag,'LineWidth',lw);
    end
    plot3(nod11(:,1),nod11(:,2),repmat(z,size(nod11,1),1),colortag,'LineWidth',lw);
    plot3(nod12(:,1),nod12(:,2),repmat(z,size(nod12,1),1),colortag,'LineWidth',lw);
    plot3(nod13(:,1),nod13(:,2),repmat(z,size(nod13,1),1),colortag,'LineWidth',lw);
    plot3(nod14(:,1),nod14(:,2),repmat(z,size(nod14,1),1),colortag,'LineWidth',lw);
    plot3(nod15(:,1),nod15(:,2),repmat(z,size(nod15,1),1),colortag,'LineWidth',lw);
end
end





