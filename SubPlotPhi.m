function [] = SubPlotPhi(mesh, params, phi, file_name, save)

if nargin < 5
    save = 0;
end

phi_final = round(phi);
ii1 = phi_final == 0;
ii2 = phi_final == 1;
ii3 = ~(ii1 | ii2);

fig1 = figure();
subplot(1,2,1)
Valve_PlotEdges(params, 1);
plot(mesh.x_mid(phi < 0.25),mesh.y_mid(phi < 0.25),'o','MarkerFaceColor','b');
hold on;
plot(mesh.x_mid(phi >= 0.25 & phi <= 0.75),mesh.y_mid(phi >= 0.25 & phi <= 0.75),'yo','MarkerFaceColor','y');
plot(mesh.x_mid(phi > 0.75),mesh.y_mid(phi > 0.75),'ro','MarkerFaceColor','r');
axis equal;

subplot(1,2,2)
Valve_PlotEdges(params, 1);
plot(mesh.x_mid(ii1),mesh.y_mid(ii1),'o','MarkerFaceColor','b');
hold on;
plot(mesh.x_mid(ii3),mesh.y_mid(ii3),'yo','MarkerFaceColor','y');
plot(mesh.x_mid(ii2),mesh.y_mid(ii2),'ro','MarkerFaceColor','r');
axis equal;

if save
    saveas(fig1, file_name);
end

end

