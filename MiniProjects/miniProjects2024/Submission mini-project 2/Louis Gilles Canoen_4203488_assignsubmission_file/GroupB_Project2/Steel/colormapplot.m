function colormapplot(nodesCoordinates,connectivity,dx_nodes,nodes_deformed_exagerated,E)
%nodesCoordinates: orignal node pos in a nb of nodes by 2 array
%connectivity table
%dx_nodes: real delta displacement of all the nodes in a nb of nodes by 2 array
%nodes_deformed_exagerated: position where the nodes are going to be
%plotted
%E: young modulus to calculate the stress



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

    total_element_length = total_element_length + elementLength;
    q_loc1 = dx_nodes(localIndices(1),:);
    q_loc2 = dx_nodes(localIndices(2),:);
    q_loc_tot = [q_loc1'; q_loc2'];

    sigma_tot(i) = E/elementLength*[-C -S C S]*q_loc_tot / 1e6 ; % stress in MPascals

end


% Plot truss with colormap
figure;
hold on;
colormap(jet); % Use jet colormap (blue to red)
% Define colormap before the loop
colorMap = jet(256);  

% Plot each bar with color based on force
for i = 1:size(connectivity, 1)
    % Get node indices
    node1 = connectivity(i,1);
    node2 = connectivity(i,2);

    % Get bar coordinates
    x = [nodes_deformed_exagerated(node1,1), nodes_deformed_exagerated(node2,1)];
    y = [nodes_deformed_exagerated(node1,2), nodes_deformed_exagerated(node2,2)];
    
    % Normalize stress values to the range [1, 256]
    sigma_min = min(sigma_tot);
    sigma_max = max(sigma_tot);
    colorIndex = round((sigma_tot(i) - sigma_min) / (sigma_max - sigma_min) * 255) + 1;
    
    % Ensure index stays within valid range
    colorIndex = max(1, min(256, colorIndex));
    
    % Assign color based on stress
    barColor = colorMap(colorIndex, :); 

    % Plot bar with corresponding color
    plot(x, y, 'LineWidth', 2.5, 'Color', barColor);
end

% Add colorbar
h = colorbar;
ylabel(h, 'Axial Stress (MPa)', 'FontSize', 16, 'FontWeight', 'bold','interpreter','Latex');
caxis([min(sigma_tot), max(sigma_tot)]);
%title('Axial Force Distribution (Compression \& Tension)','interpreter','Latex');
xlabel('X Position','interpreter','Latex');
ylabel('Y Position','interpreter','Latex');
grid on;
axis equal;

% Get the coordinates of node 5
node5_x = nodes_deformed_exagerated(5,1);
node5_y = nodes_deformed_exagerated(5,2);

% Define arrow properties
arrowLength = 2; % Adjust length as needed
arrow_dx = 0;  % No horizontal movement
arrow_dy = -arrowLength;  % Downward movement

% Plot the arrow
quiver(node5_x, 1.6*node5_y, arrow_dx, arrow_dy, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);
fontsize(16,"points")
end