% Code for only one pair of (n,h) input. See "Aluminum_alloy_2" for both sweeps through n and h.

clear all
close all
clc

% Physical properties (Aluminum alloy 6061)
E = 69e9;
A = 70e-4; % Max cross-sectional area chosen
rho = 2700;

% Geometry (best parameters for highest first frequency)
h = 5.8; %5.8
n = 1.06; %1.06

multi = 1000; % multiplication factor for displacement plot

% Applied force
f = 50e3;
f_vect = zeros(14, 1);
f_vect(2) = -f;

% functions for node height
best = @(x, y, n) y*(1 - (x-10)^2/121)^n;

% Test for nodes 10 and 11 (distance = 5/11 in the study) --> This was for some testing
distance = 10/21; % distance = 0 means nodes 9 and 10 are on node 5, distance = 1 means they are on top of nodes 1 and 3.

% Geometry
NodesCoord = [0, 0;    %1
              -1, 0;    %2
              20, 0;   %3
              21, 0;   %4
              10, h;  %5
              20/3, 0;  %6
              40/3, 0;  %7
              -0.5, best(-0.5,h,n);  %8
              (1-distance)*10 + distance*(-0.5), best((1-distance)*10 + distance*(-0.5),h,n);  %9 (5, ...)
              (1-distance)*10 + distance*(20.5), best((1-distance)*10 + distance*(20.5),h,n); %10 (15, ...)
              20.5, best(20.5,h,n)]; %11

connectivity = [1, 2;  %1
                3, 4;  %2
                1, 6;  %3
                6, 7;  %4
                3, 7;  %5
                2, 8; %6
                1, 8; %7
                1, 9; %8
                6, 9; %9
                5, 6; %10
                5, 7; %11
                7, 10; % 12
                3, 10; %13
                3, 11; %14
                4, 11; %15
                8, 9; %16
                5, 9; %17
                5, 10; %18
                10, 11]; %19

% Volume used
Volume = A * sum(vecnorm(NodesCoord(connectivity(:,1), :) - NodesCoord(connectivity(:,2), :), 2, 2));

numberOfNodes = size(NodesCoord,1);
GDof = 2*numberOfNodes;

K = formStiffness2Dtruss(GDof, connectivity, NodesCoord, E*A);
M = formConsistentMass2Dtruss(GDof, connectivity, NodesCoord, rho*A); % We only used the consistent mass matrix for this case

%% Boundary condition + solving the eigenvalue problem

FixedDOF = [1, 2, 3, 4, 5, 6, 7, 8]';

ActiveDOF = setdiff((1:GDof)', FixedDOF);
Kred = K(ActiveDOF, ActiveDOF);
Mred = M(ActiveDOF, ActiveDOF);

[modal_matrix, omega_method1] = computeFrequenciesAndModes(GDof, FixedDOF, K, M, 0);

value = omega_method1(1)/(2*pi);

modes_normalized = zeros(size(modal_matrix));

for i = 1:size(modal_matrix, 2)
    norm = sqrt(modal_matrix(:, i))' * Mred * modal_matrix(:, i);
    modes_normalized(:, i) = modal_matrix(:, i)./norm;
end

% Extract the first 5 natural frequencies (in Hz)
freq_5_method1 = omega_method1(1:length(omega_method1))/(2*pi);
FrequencyTable = table((1:length(omega_method1))', freq_5_method1, ...
    'VariableNames', {'ModeNumber', 'Consistent_Hz'});
% Display the table
disp(FrequencyTable);

%% Static

% Finding the delta x displacements
delta_x = Kred\f_vect;

% Reshaping the 14x1 vector to a 7x2 matrix 
delta_x_corr = reshape(delta_x, 2, [])';

% Expand delta_x to include zero displacements for fixed nodes
full_delta_x = zeros(numberOfNodes, 2);  % Initialize an 8x2 matrix

% Assign the computed displacements to nodes 5-8
full_delta_x(5:end, :) = delta_x_corr;  

% Create the deformed coordinate matrix
displacement = NodesCoord + full_delta_x;
displacement_plot = NodesCoord + multi*full_delta_x; % exagerated displacements for plot later on

%% Compression/Tension Colormap

nodeCoordinateX = NodesCoord(:,1);
nodeCoordinateY = NodesCoord(:,2);
total_element_length = 0;
sigma_tot = zeros(1,size(connectivity,1));

for i = 1:size(connectivity, 1)
    localIndices = connectivity(i,:);   
    xa = nodeCoordinateX(localIndices(2))-nodeCoordinateX(localIndices(1));
    ya = nodeCoordinateY(localIndices(2))-nodeCoordinateY(localIndices(1));
    elementLength = sqrt(xa*xa+ya*ya);
    C = xa/elementLength;
    S = ya/elementLength;

    total_element_length = total_element_length + elementLength;
    q_loc1 = full_delta_x(localIndices(1),:);
    q_loc2 = full_delta_x(localIndices(2),:);
    q_loc_tot = [q_loc1'; q_loc2'];

    sigma_tot(i) = E/elementLength*[-C -S C S]*q_loc_tot / 1e6 ; % stress in MPa
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
    x = [displacement_plot(node1,1), displacement_plot(node2,1)];
    y = [displacement_plot(node1,2), displacement_plot(node2,2)];
    
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
ylabel(h, 'Axial Stress [MPa]', 'FontSize', 12, ...
    'FontWeight', 'bold', 'Interpreter', 'latex');
clim([min(sigma_tot), max(sigma_tot)]);

%title('Axial Force Distribution (Compression \& Tension)', 'Interpreter', 'latex');
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
xlim([-1 21])
ylim([-4 10])

% Get the coordinates of node 5
node5_x = displacement_plot(5,1);
node5_y = displacement_plot(5,2);

% Define arrow properties (for visualization of the applied point force)
arrowLength = 2; % Adjust length as needed
arrow_dx = 0;  % No horizontal movement
arrow_dy = -arrowLength;  % Downward movement

% Plot the arrow
quiver(node5_x, 1.5*node5_y, arrow_dx, arrow_dy, ...
    'k', 'LineWidth', 2, 'MaxHeadSize', 1);

% Use colormap colors at the edges (min = index 1, max = index 256)
minColor = colorMap(1, :);
maxColor = colorMap(end, :);

% Create dummy plots with those colors
dummyMax = plot(NaN, NaN, '-', 'LineWidth', 2.5, 'Color', maxColor);
dummyMin = plot(NaN, NaN, '-', 'LineWidth', 2.5, 'Color', minColor);

% Create legend with stress info
legendStr = {
    sprintf('Max stress: %.2f MPa', max(sigma_tot)), ...
    sprintf('Min stress: %.2f MPa', min(sigma_tot))
};

legend([dummyMax, dummyMin], legendStr, ...
    'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 12);

% Printing important information
fprintf('f_1 = %.3f Hz, max stress = %.3f MPa, min_stress = %.3f MPa, volume = %.3f m^3, max displacement %.3f mm \n', value, max(sigma_tot), min(sigma_tot), Volume, 1000*max(vecnorm(full_delta_x,2,2)));

%% perfect plot (2D truss bridge plot)

figure
hold on
% Plot each bar with color based on force
for i = 1:size(connectivity, 1)
    % Get node indices
    node1 = connectivity(i,1);
    node2 = connectivity(i,2);

    % Get bar coordinates
    x = [displacement_plot(node1,1), displacement_plot(node2,1)];
    y = [displacement_plot(node1,2), displacement_plot(node2,2)];

    % Plot bar with corresponding color
    plot(x, y, 'LineWidth', 2.5, 'Color', [0.4, 0.4, 0.4]);
end

for i = 1:size(NodesCoord, 1)

    % Get bar coordinates
    x = NodesCoord(i,1);
    y = NodesCoord(i,2);

    % Plot bar with corresponding color
    plot(x, y, 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1]);
end

% dummy plots for legend
p1 = plot(NaN, NaN, 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1]);
p2 = plot(NaN, NaN, '-', 'LineWidth', 2.5, 'Color', [0.4, 0.4, 0.4]);
xlim([-1, 21])
ylim([-5 10])
legend([p1, p2], {'Nodes', 'Bars'}, 'Interpreter', 'latex', 'Fontsize', 13)
axis off

%% Modal Analysis

modeNumber = 3; % <--------------------- CHOOSE HERE THE MODE SHAPE TO DISPLAY

us = 1:2:2*numberOfNodes-1;
vs = 2:2:2*numberOfNodes;
full_modes = zeros(GDof,length(omega_method1));
activeDof = setdiff(transpose((1:GDof)), FixedDOF);
full_modes(activeDof,:) = modal_matrix;
XX = full_modes(us,modeNumber); 
YY = full_modes(vs,modeNumber);
dispNorm = max(sqrt(XX.^2+YY.^2));

scaleFact = 300*dispNorm; % <------------------ CHOOSE HERE THE SCALING FACTOR

figure
set(gca, 'Position', [0.1, 0.15, 0.85, 0.6]);
hold on

% Plot phantom (undeformed) structure in light gray
for i = 1:size(connectivity, 1)
    node1 = connectivity(i,1);
    node2 = connectivity(i,2);

    x = [NodesCoord(node1,1), NodesCoord(node2,1)];
    y = [NodesCoord(node1,2), NodesCoord(node2,2)];

    plot(x, y, '-', 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5, 0.5]);
end

% Prepare arrays for deformed shapes (positive and negative)
X_pos = NodesCoord(:,1) + scaleFact * XX;
Y_pos = NodesCoord(:,2) + scaleFact * YY;

X_neg = NodesCoord(:,1) - scaleFact * XX;
Y_neg = NodesCoord(:,2) - scaleFact * YY;

% Fill area between positive and negative deformed shapes
for i = 1:size(connectivity, 1)
    node1 = connectivity(i,1);
    node2 = connectivity(i,2);

    X_fill = [X_pos(node1), X_pos(node2), X_neg(node2), X_neg(node1)];
    Y_fill = [Y_pos(node1), Y_pos(node2), Y_neg(node2), Y_neg(node1)];

    fill(X_fill, Y_fill, [0.8, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

% Plot positive deformed shape in red
for i = 1:size(connectivity, 1)
    node1 = connectivity(i,1);
    node2 = connectivity(i,2);

    x = [X_pos(node1), X_pos(node2)];
    y = [Y_pos(node1), Y_pos(node2)];

    plot(x, y, 'r-', 'LineWidth', 2);
end

% Plot negative deformed shape in blue
for i = 1:size(connectivity, 1)
    node1 = connectivity(i,1);
    node2 = connectivity(i,2);

    x = [X_neg(node1), X_neg(node2)];
    y = [Y_neg(node1), Y_neg(node2)];

    plot(x, y, 'b-', 'LineWidth', 2);
end

% Create dummy handles for legend
h1 = plot(NaN, NaN, '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5 0.5]);
h2 = fill(NaN, NaN, [0.8 0.3 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h3 = plot(NaN, NaN, 'r-', 'LineWidth', 2.5);
h4 = plot(NaN, NaN, 'b-', 'LineWidth', 2.5);

axis equal
title(sprintf('Mode %d Shape (%.2f Hz)', modeNumber, freq_5_method1(modeNumber)), ...
      'Interpreter', 'latex', 'FontSize', 14)
xlabel('$x$ [m]', 'Interpreter', 'latex')
ylabel('$y$ [m]', 'Interpreter', 'latex')
ylim([-6 8])
axis off

%% MOVIE
% Uncomment from here down to plot a small 5 second movie of the mode shape moving!
% 
% % Animate vibration of the truss for selected mode
% figure;
% axis equal
% hold on
% xlim([-1 21])
% ylim([-5 10])
% title(sprintf('Mode %d vibration shape', modeNumber), 'Interpreter', 'latex', 'FontSize', 14);
% xlabel('x [m]', 'Interpreter', 'latex');
% ylabel('y [m]', 'Interpreter', 'latex');
% 
% % Plot phantom structure (undeformed, gray)
% for i = 1:size(connectivity, 1)
%     node1 = connectivity(i,1);
%     node2 = connectivity(i,2);
%     x = [NodesCoord(node1,1), NodesCoord(node2,1)];
%     y = [NodesCoord(node1,2), NodesCoord(node2,2)];
%     plot(x, y, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); 
% end
% 
% 
% 
% % Set up animated plot (initial frame)
% lineHandles = gobjects(size(connectivity,1),1);
% for i = 1:size(connectivity, 1)
%     lineHandles(i) = plot(NaN, NaN, 'r-', 'LineWidth', 2.5);
% end
% 
% % Vibration animation
% fps = 60;
% duration = 5; % seconds
% t = linspace(0, duration, fps*duration);
% 
% % FOR VIDEO
% % % Setup video writer
% % videoName = sprintf('mode_%d_vibration.mp4', modeNumber);
% % v = VideoWriter(videoName, 'MPEG-4');
% % v.FrameRate = fps;
% % open(v);
% % FOR VIDEO
% 
% for i = 1:length(t)
%     phase = sin(2*pi*t(i));
%     Xdisp = NodesCoord(:,1) + scaleFact * XX * phase;
%     Ydisp = NodesCoord(:,2) + scaleFact * YY * phase;
% 
%     for j = 1:size(connectivity, 1)
%         n1 = connectivity(j,1);
%         n2 = connectivity(j,2);
%         set(lineHandles(j), 'XData', [Xdisp(n1), Xdisp(n2)], ...
%                             'YData', [Ydisp(n1), Ydisp(n2)]);
%     end
% 
%     drawnow;
%     % writeVideo(v, getframe(gcf));  % Save current frame % NEED FOR VIDEO
% end
% 
% % close(v); % NEED FOR VIDEO
% % fprintf('✅ Video saved as "%s"\n', videoName); % NEED FOR VIDEO

% figure
% set(gca, 'Position', [0, 0, 1, 1]); % fill figure window
% 
% % Create dummy handles
% h1 = plot(NaN, NaN, '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5 0.5]); hold on
% h2 = fill(NaN, NaN, [0.8 0.3 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% h3 = plot(NaN, NaN, 'r-', 'LineWidth', 2.5);
% h4 = plot(NaN, NaN, 'b-', 'LineWidth', 2.5);
% 
% % Create legend
% lgd = legend([h1 h2 h3 h4], {'Undeformed', 'Sweep', 'Mode +', 'Mode -'}, ...
%        'Interpreter', 'latex', 'FontSize', 11, ...
%        'NumColumns', 1, 'Orientation', 'horizontal', 'Fontsize', 30);
% 
% % Position legend manually in the center
% set(lgd, 'Position', [0, 0, 1, 1]); % [left, bottom, width, height]
% 
% axis off