function [] = Valve_PlotEdges(z)

load(fullfile('Valve_Data', 'param'));

nod1 = [
    0, -y_air/2;      % Air region
    0, y_air/2;
    x_air, y_air/2;
    x_air, -y_air/2;
    0, -y_air/2; 
    ] ;

nod2 = [
    x_channel + x1 + 2*t, (y3/2) + t + y2;        % Conductor 1 (upper)
    x_channel + x1 + 2*t, (y3/2) + t;
    x_channel + x1 + 2*t + x2, (y3/2) + t;
    x_channel + x1 + 2*t + x2, (y3/2) + t + y2;
    x_channel + x1 + 2*t, (y3/2) + t + y2;
    ] ;


nod3 = [
    (x_channel + x1 +2*t), -((y3/2) + t + y2);        % Conductor 2 (lower)
    (x_channel + x1 +2*t), -((y3/2) + t);
    (x_channel + x1 +2*t + x2), -((y3/2) + t);
    (x_channel + x1 +2*t + x2), -((y3/2) + t + y2);
    (x_channel + x1 +2*t), -((y3/2) + t + y2); 
    ] ;

nod4 = [
    (x_channel),  -(y3/2 + 2*t + y2 - y4 - t);        % Plunger
    (x_channel + x1),  -(y3/2 + 2*t + y2 - y4 - t);
    (x_channel + x1),  -(y3/2 + 2*t + y2 - y4 - t - y_piston);
    (x_channel),  -(y3/2 + 2*t + y2 - y4 - t - y_piston);
    (x_channel),  -(y3/2 + 2*t + y2 - y4 - t);
    ] ;

nod5 = [
    x_channel, (y3/2) + y2 + 2*t - y4;        % Iron
    x_channel, (y3/2) + y2 + 2*t + y1;
    x_channel + x1 + x4 + w_m + t + x3, (y3/2) + y2 + 2*t + y1;
    x_channel + x1 + x4 + w_m + t + x3, -((y3/2) + y2 + 2*t + y1);
    x_channel, -((y3/2) + y2 + 2*t + y1);
    x_channel, -((y3/2) + y2 + 2*t - y4);
    x_channel + x1, -((y3/2) + y2 + 2*t - y4);
    x_channel + x1, -((y3/2) + y2 + 2*t);
    x_channel + x1 + x4 + w_m + 2*t , -((y3/2) + y2 + 2*t);
    x_channel + x1 + x4 + w_m + 2*t , -((y3/2));
    x_channel + x1 + x4 + w_m - x2 , -((y3/2));
    x_channel + x1 + x4 + w_m - x2 , (y3/2);
    x_channel + x1 + x4 + w_m + 2*t , (y3/2);
    x_channel + x1 + x4 + w_m + 2*t , (y3/2) + 2*t + y2;
    x_channel + x1 , (y3/2) +  2*t + y2;
    x_channel + x1 , (y3/2) + 2*t + y2 - y4;
    x_channel, (y3/2) + y2 + 2*t - y4;   
    ] ;

if nargin < 1
    hold on;
    plot(nod1(:,1),nod1(:,2),'k-');    
    plot(nod2(:,1),nod2(:,2),'k-');
    plot(nod3(:,1),nod3(:,2),'k-');
    plot(nod4(:,1),nod4(:,2),'k-');
    plot(nod5(:,1),nod5(:,2),'k-');
    axis equal
else
    hold on;
    plot3(nod1(:,1),nod1(:,2),repmat(z,size(nod1,1),1),'k-');    
    plot3(nod2(:,1),nod2(:,2),repmat(z,size(nod2,1),1),'k-');
    plot3(nod3(:,1),nod3(:,2),repmat(z,size(nod3,1),1),'k-');
    plot3(nod4(:,1),nod4(:,2),repmat(z,size(nod4,1),1),'k-');
    plot3(nod5(:,1),nod5(:,2),repmat(z,size(nod5,1),1),'k-');
end
end





