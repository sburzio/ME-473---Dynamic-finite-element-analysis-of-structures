function [nodesCoordinates,connectivity] = Howe(complexity,base)
%Creates a howe style truss and gives you node coordinates and connectivity
%table
%   complexity = what order of truss do you want (>= 0, integer values)
%   a = height of max point
%   base = base node coordinates (1 - 5) with node 5 being the peak, 1 and
%   2 being the leftmost base points, and 3 and 4 being rightmost base
%   points (points 1 and 3 are inner, 2 and 4 are outer)

if complexity == 0
    center_x = mean([base(1,1),base(3,1)]); % center x point
    nodesCoordinates = [base;center_x,0];
    connectivity = [1 2; 3 4; 1 6; 6 3; 6 5; 2 5; 4 5];
else
    num_nodes = 8 + 4*(complexity - 1); % total number of nodes
    num_nodes_bottom = 5 + 2*(complexity - 1); % number of nodes on the bottom of the truss
    
    center_x = mean([base(1,1),base(3,1)]); % center x point
    
    % create vector of nodes along the bottom of the truss (not including base
    % points)
    nodes_bot_inner = zeros(num_nodes_bottom-2,2);
    nodes_bot_inner(:,1) = linspace(base(1,1),base(3,1),complexity*2+1);
    nodes_bot_inner = nodes_bot_inner(2:end-1,:);
    
    % find the angle between base points and peak, as well as distance
    ang = atan(base(5,2)/(center_x-base(2,1)));
    dist_top_tot = norm(base(5,:)-base(2,:));
    dist_top = dist_top_tot./(complexity+1);
    
    % Find the nodes on the outer right and left edges
    nodes_top_left = zeros(complexity,2);
    nodes_top_right = zeros(complexity,2);
    for i = 1:complexity
        if i == 1
            nodes_top_left(i,:) = [base(1,1), base(2,2)+ tan(ang).*(base(1,1)-base(2,1))]; 
            nodes_top_right(i,:) = [base(3,1), base(4,2)+ tan(ang).*(base(4,1)-base(3,1))]; 
        else
            nodes_top_left(i,:) = [nodes_bot_inner(i-1,1), base(2,2)+ tan(ang).*(nodes_bot_inner(i-1,1)-base(2,1))];
            nodes_top_right(i,:) = [nodes_bot_inner(end-(i-2),1), base(4,2)+ tan(ang).*(base(4,1)-nodes_bot_inner(end-(i-2),1))];
        end
    end
    
    nodes_top_right = nodes_top_right(end:-1:1,:);
    
    % create overall vector of nodes using the 5 base nodes, followed by the
    % inner bottom nodes, followed by top left edge nodes, followed by top
    % right edge nodes
    nodesCoordinates = [base; nodes_bot_inner; nodes_top_left; nodes_top_right]; 
    
    % create connectivity table of base
    connect_bot_base = [1 2; 3 4];
    
    % create connectivity table for inner bottom beams
    connect_bot_inner = zeros(2*complexity,2);
    for i = 1:2*complexity
        if i == 1
            connect_bot_inner(i,:) = [1,6];
        elseif i == 2*complexity
            connect_bot_inner(i,:) = [6+2*(complexity-1),3];
        else
            connect_bot_inner(i,:) = [old,old+1];
        end
        old = connect_bot_inner(i,2);
    end
    
    last_bot = connect_bot_inner(i-1,2);
    
    % create connectivity table for upper left beams
    connect_left = zeros(1+complexity,2);
    for i = 1:1+complexity
        if i == 1
            connect_left(i,:) = [2,last_bot+1];
        elseif i == 1+complexity
            connect_left(i,:) = [last_bot+complexity,5];
        else
            connect_left(i,:) = [old,old+1];
        end
        old = connect_left(i,2);
    end
    
    last_left = connect_left(i-1,2);
    
    
    % create connectivity table for upper right beams
    connect_right = zeros(1+complexity,2);
    for i = 1:1+complexity
        if i == 1
            connect_right(i,:) = [5,last_left+1];
        elseif i == 1+complexity
            connect_right(i,:) = [last_left+complexity,4];
        else
            connect_right(i,:) = [old,old+1];
        end
        old = connect_right(i,2);
    end
    
    last_right = connect_right(i-1,2);
    
    % create connectivity table for inner beams, left side
    connect_inner_left = zeros(2*complexity,2);
    bot_curr = 6;
    top_curr = last_bot+1;
    vert = true; % always start with vertical bar
    for i = 1:2*complexity
        if i == 1 
            connect_inner_left(i,:) = [1,top_curr];
            top_curr = top_curr+1;
        else
            if vert == true % if the last bar was vertical, the next one is at an angle
                connect_inner_left(i,:) = [old,bot_curr];
                bot_curr = bot_curr +1;
                vert = false;
            else
                connect_inner_left(i,:) = [old,top_curr];
                top_curr = top_curr +1;
                vert = true;
            end
        end
        old = connect_inner_left(i,2);
    end
    
    % create connectivity table for inner beams, right side
    connect_inner_right = zeros(2*complexity,2);
    bot_curr = connect_inner_left(end,2) + 1;
    top_curr = last_left + 1;
    vert = true; % always start with vertical bar
    for i = 1:2*complexity+1 % include the middle bar in this group
        if i == 1 
            connect_inner_right(i,:) = [5,connect_inner_left(end,2)];
        elseif i == 2*complexity+1
            connect_inner_right(i,:) = [connect_inner_right(i-1,2),3];
        else
            if vert == 1 % if the last bar was vertical, the next one is at an angle
                connect_inner_right(i,:) = [old,top_curr];
                top_curr = top_curr +1;
                vert = false;
            else
                connect_inner_right(i,:) = [old,bot_curr];
                bot_curr = bot_curr +1;
                vert = true;
            end
        end
        old = connect_inner_right(i,2);
    end
    
    connectivity = [connect_bot_base; connect_bot_inner; connect_left; connect_right; connect_inner_left; connect_inner_right];

end

end