function draw2Dtruss(nodesCoordinates, connectivity, varargin)
% draw2Dtruss Draws a 2D truss structure with node and element labels
% 
% Usage:
%   draw2Dtruss(nodesCoordinates, connectivity)
%   draw2Dtruss(nodesCoordinates, connectivity, 'LineColor', 'r', 'Title', 'My Truss')
%
% Inputs:
%   nodesCoordinates - Nx2 matrix of node coordinates [x, y]
%   connectivity     - Mx2 matrix of node indices defining truss elements
%   varargin         - (Optional) name-value pairs:
%                      'LineColor' (default: 'b')
%                      'LineWidth' (default: 1)
%                      'MarkerSize' (default: 8)
%                      'ShowNodeNumbers' (default: false)
%                      'ShowBeamNumbers' (default: false)
%                      'Title' (default: '2D Truss Structure')

    % Parse input options
    p = inputParser;
    addParameter(p, 'LineColor', 'b');
    addParameter(p, 'LineWidth', 1);
    addParameter(p, 'MarkerSize', 15);
    addParameter(p, 'ShowNodeNumbers', false);
    addParameter(p, 'ShowBeamNumbers', false);
    addParameter(p, 'Title', '2D Truss Structure');
    addParameter(p, 'FontSize', 14)
    parse(p, varargin{:});
    
    % Extract parsed options
    lineColor = p.Results.LineColor;
    lineWidth = p.Results.LineWidth;
    markerSize = p.Results.MarkerSize;
    showNodeNumbers = p.Results.ShowNodeNumbers;
    showBeamNumbers = p.Results.ShowBeamNumbers;
    fontSize = p.Results.FontSize;
    titleText = p.Results.Title;

    % figure; 
    hold on; axis equal;
    xlabel('X'); ylabel('Y');
    title(titleText, 'Interpreter', 'none');

    % Plot truss elements
    for i = 1:size(connectivity, 1)
        n1 = connectivity(i, 1);
        n2 = connectivity(i, 2);
        x = [nodesCoordinates(n1, 1), nodesCoordinates(n2, 1)];
        y = [nodesCoordinates(n1, 2), nodesCoordinates(n2, 2)];
        plot(x, y, '-', 'Color', lineColor, 'LineWidth', lineWidth);

        % Midpoint of element
        midX = mean(x);
        midY = mean(y);

        % % Plot white triangle with black edge at midpoint
        % plot(midX, midY, '^', ...
        %     'MarkerSize', markerSize * 0.8, ...
        %     'MarkerEdgeColor', 'k', ...
        %     'MarkerFaceColor', 'w');

        % Label the truss element number (larger, black font)
        if showBeamNumbers
            text(midX, midY, sprintf('%d', i), ...
                'Color', 'b', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end

    % Plot nodes with white interior circles
    plot(nodesCoordinates(:,1), nodesCoordinates(:,2), '.', ...
         'MarkerSize', markerSize, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

    % Node numbers inside the circle
    if showNodeNumbers
        for i = 1:size(nodesCoordinates, 1)
            text(nodesCoordinates(i,1), nodesCoordinates(i,2), ...
                sprintf('%d', i), ...
                'FontSize', fontSize, 'Color', 'k', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top');
        end
    end

    grid on;
end