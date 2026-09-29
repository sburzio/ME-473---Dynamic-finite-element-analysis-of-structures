function [handle_leg] = noNumPlot(nodesCoordinates,connectivity,color,fignum)
% nodesCoordinates: orignal node pos in a nb of nodes by 2 array
% connectivity table


nodeCoordinateX = nodesCoordinates(:,1);
nodeCoordinateY = nodesCoordinates(:,2);
total_element_length = 0;

for i = 1:size(connectivity, 1)
    localIndices = connectivity(i,:);   
    xa = nodeCoordinateX(localIndices(2))-nodeCoordinateX(localIndices(1));
    ya = nodeCoordinateY(localIndices(2))-nodeCoordinateY(localIndices(1));
    elementLength = sqrt(xa*xa+ya*ya);
    C = xa/elementLength;
    S = ya/elementLength;

end

% Plot truss
figure(fignum);
hold on;

% Plot each bar 
for i = 1:size(connectivity, 1)
    % Get node indices
    node1 = connectivity(i,1);
    node2 = connectivity(i,2);

    % Get bar coordinates
    x = [nodesCoordinates(node1,1), nodesCoordinates(node2,1)];
    y = [nodesCoordinates(node1,2), nodesCoordinates(node2,2)];

    % Plot bar with corresponding color
    handle_leg = plot(x, y, 'LineWidth', 2.5, 'Color', color);

end

% Plot nodes
scatter(nodesCoordinates(:,1),nodesCoordinates(:,2),[],color,'filled','o') 

hold off
xlabel('X Position','interpreter','latex');
ylabel('Y Position','interpreter','latex');
axis equal;
fontsize(16,"points")
%axis off
grid on
grid minor


end