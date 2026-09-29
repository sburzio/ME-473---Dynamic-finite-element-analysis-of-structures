clear
close all
clc

%change the name in case its different
load('homework2_data_2.mat')

%plot variable
nb_node_inter=[4 6 8 10 14 20 40]
a=linspace(1,20,nb_a)

%surface plot (didnt used it in the report since the 4-su8blot was enough
%however it is pratical to check the data
[X,Y] = meshgrid(a,nb_node_inter);
figure
surf(X,Y,maxdisplacement')
title('max displacement')
xlabel('height of the node 5 (a)')
ylabel('number of nodes of intermediate nodes on each side')
figure
surf(X,Y,maxabsstress')
title('max absolute stress')
xlabel('height of the node 5 (a)')
ylabel('number of nodes of intermediate nodes on each side')
figure
surf(X,Y,volume')
title('volume')
xlabel('height of the node 5 (a)')
ylabel('number of nodes of intermediate nodes on each side')
figure
surf(X,Y,frequency1')
title('natural frequency')
xlabel('height of the node 5 (a)')
ylabel('number of nodes of intermediate nodes on each side')


%4-subplot plot
figure('Position',[100 100 900 400]);
subplot(2,2,1)
plot(a,volume)
ylabel('\textbf{Total volume of the beams}')
xlabel('\textbf{a}')
xlim([1 20])
xline(2.5,'Label','a=2.5m','Interpreter','latex')
subplot(2,2,2)
plot(a,frequency1)
ylabel('\textbf{First frequency}')
xlabel('\textbf{a}')
xlim([1 20])
xline(2.5,'Label','a=2.5m','Interpreter','latex')
subplot(2,2,3)
plot(a,maxdisplacement)
ylabel('\textbf{Max displacement (mm)}')
xlabel('\textbf{a}')
xlim([1 20])
xline(2.5,'Label','a=2.5m','Interpreter','latex')
subplot(2,2,4)
plot(a,maxabsstress)
xlabel('\textbf{a}')
ylabel('\textbf{Max absolute stress}')
xlim([1 20])
xline(2.5,'Label','a=2.5m','Interpreter','latex')
fontsize(scale=1.2)
legend(['nb inter. node: '+string(nb_node_inter)],'Location','best')
exportgraphics(gcf,'wood_results.pdf','ContentType','vector')