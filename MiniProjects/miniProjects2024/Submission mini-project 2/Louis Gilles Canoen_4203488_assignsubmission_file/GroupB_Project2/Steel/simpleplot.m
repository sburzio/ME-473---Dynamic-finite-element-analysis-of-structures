function simpleplot(nodesCoordinates,connectivity,complexity,buffer,bar_plot)
% nodesCoordinates: orignal node pos in a nb of nodes by 2 array
% connectivity table
% buffer = desired spacing between nodes and their labels
% bar_plot = do you also want to plot numbers of bars? true = yes

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

centers = zeros(size(connectivity,1),2); % empty vector to hold center position of each bar
labels_bars = cell(size(connectivity,1),1);
ct_bar = 1;

for i = 1:size(connectivity,1)
    labels_bars{i} = num2str(ct_bar);
    ct_bar = ct_bar+1;
end

% Plot truss
figure;
hold on;

% Plot each bar with color based on force
for i = 1:size(connectivity, 1)
    % Get node indices
    node1 = connectivity(i,1);
    node2 = connectivity(i,2);

    % Get bar coordinates
    x = [nodesCoordinates(node1,1), nodesCoordinates(node2,1)];
    y = [nodesCoordinates(node1,2), nodesCoordinates(node2,2)];

    % Plot bar with corresponding color
    plot(x, y, 'LineWidth', 2.5, 'Color', [79/256,79/256,67/256]);

    % get out center point
    centers(i,:) = [mean(x),mean(y)];
    if bar_plot == true
        scatter(centers(i,1), centers(i,2),150,'w','filled','o') 
    end
end

labels = cell(length(nodesCoordinates),1);
for i = 1:length(nodesCoordinates)
    labels{i} = num2str(i);
end

% label nodes
if complexity ~= 0
    bot_end = 5+1+2*(complexity-1);
else
    bot_end = 6;
end
top_left_end = bot_end+complexity;
scatter(nodesCoordinates(:,1),nodesCoordinates(:,2),[],'b','filled','o') 
labelpoints(nodesCoordinates(1:4,1),nodesCoordinates(1:4,2),labels(1:4),'S',buffer,1,'interpreter','Latex','Color','b')
labelpoints(nodesCoordinates(5,1),nodesCoordinates(5,2),labels(5),'N',buffer,1,'interpreter','Latex','Color','b')
labelpoints(nodesCoordinates((5+1):bot_end,1),nodesCoordinates((5+1):bot_end,2),labels((5+1):bot_end),'S',buffer,1,'interpreter','Latex','Color','b')
labelpoints(nodesCoordinates(bot_end+1:bot_end+complexity,1),nodesCoordinates(bot_end+1:bot_end+complexity,2),labels(bot_end+1:bot_end+complexity),'NW',buffer,1,'interpreter','Latex','Color','b')
labelpoints(nodesCoordinates(top_left_end+1:top_left_end+complexity,1),nodesCoordinates(top_left_end+1:top_left_end+complexity,2),labels(top_left_end+1:top_left_end+complexity),'NE',buffer,1,'interpreter','Latex','Color','b')

% label bars
if bar_plot == true
    labelpoints(centers(:,1),centers(:,2),labels_bars,'center',buffer,1,'interpreter','Latex','Color',[79/256,79/256,67/256])
end

hold off
%xlabel('X Position');
%ylabel('Y Position');
axis equal;
fontsize(16,"points")
axis off

% Get the coordinates of node 5
node5_x = nodesCoordinates(5,1);
node5_y = nodesCoordinates(5,2);

end