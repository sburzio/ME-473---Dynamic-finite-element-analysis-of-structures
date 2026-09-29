% Code for the study of n and h on the bridge's performance

clear all
close all
clc

% Physical properties (aluminum alloy 6061)
E = 69e9;
A = 70e-4;
rho = 2700;

% Applied force
f = 50e3;
f_vect = zeros(14, 1);
f_vect(2) = -f;

% Initialisation
VolumeN_H = zeros(201, 100);
maxStressN_H = zeros(201, 100);
minStressN_H = zeros(201, 100);
freqN_H = zeros(201, 100);
sigma_tot = zeros(1, 19);
n_index = 0;

% function for node height
best = @(x, y, n) y*(1 - (x-10)^2/121)^n;

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

% Test for nodes 10 and 11 (distance = 5/11 in the study) --> This was for some testing
distance = 10/21; % distance = 0 means nodes 9 and 10 are on node 5, distance = 1 means they are on top of nodes 1 and 3.

n_vect = 0:0.01:2;
h_vect = 0.1:0.1:10; % h can't start at 0 or else K will become singular!

for n = n_vect
h_index = 0;
n_index = n_index + 1;
n % Dispaly n to see progress of the for-loop

for h = h_vect
h_index = h_index + 1;

% Geometry (NodesCoord is defined inside the loop as coordinates change all the time)
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

numberOfNodes = size(NodesCoord,1);
GDof = 2*numberOfNodes;

K = formStiffness2Dtruss(GDof, connectivity, NodesCoord, E*A);
M = formConsistentMass2Dtruss(GDof, connectivity, NodesCoord, rho*A);

%% Boundary condition + solving the eigenvalue problem

FixedDOF = [1, 2, 3, 4, 5, 6, 7, 8]';

ActiveDOF = setdiff((1:GDof)', FixedDOF);
Kred = K(ActiveDOF, ActiveDOF);
Mred = M(ActiveDOF, ActiveDOF);

[modal_matrix, omega_method1] = computeFrequenciesAndModes(GDof, FixedDOF, K, M, 0);

value = omega_method1(1)/(2*pi);

% Frequency
freqN_H(n_index, h_index) = value; % Getting the first natural frequency of the (n,h) pair

% Volume
Volume = A * sum(vecnorm(NodesCoord(connectivity(:,1), :) - NodesCoord(connectivity(:,2), :), 2, 2));
VolumeN_H(n_index, h_index) = Volume; % Getting the volume of the (n,h) pair

% Stresses

% Expand delta_x to include zero displacements for fixed nodes
full_delta_x = zeros(numberOfNodes, 2);  % Initialize an 8x2 matrix

delta_x = Kred\f_vect;

delta_x_corr = reshape(delta_x, 2, [])';

% Assign the computed displacements to nodes 5-8
full_delta_x(5:end, :) = delta_x_corr; 

nodeCoordinateX = NodesCoord(:,1);
nodeCoordinateY = NodesCoord(:,2);
total_element_length = 0;

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

    sigma_tot(i) = E/elementLength*[-C -S C S]*q_loc_tot / 1e6 ; % stress in MPascals
end

MaxStress = max(sigma_tot);
maxStressN_H(n_index, h_index) = MaxStress; % Getting the max stress of the (n,h) pair

MinStress = min(sigma_tot);
minStressN_H(n_index, h_index) = MinStress; % Getting the min stress of the (n,h) pair

end % end of h-loop

end % end of n-loop
% ----------------------------------

%% PLOTS

% Create a mesh grid for n and h
[H, N] = meshgrid(h_vect, n_vect);  % Note: h corresponds to rows, n corresponds to columns

% Get max freq with respect to h and n
[max_freq_val, freq_row_idx] = max(freqN_H(:)); % Get the max value and its linear index
[row, col] = ind2sub(size(freqN_H), freq_row_idx); % Convert linear index to row/col



% Getting the best value of frequency and the indices of n and h
dummy_n = (0:0.01:2);
dummy_h = (0.1:0.1:10);

best_n = dummy_n(row);
best_h = dummy_h(col);
best_freq = max_freq_val;

fprintf('Best n = %.3f, best h = %.3f, best freq = %.3f \n', best_n, best_h, best_freq);

threshold = 50; % Bonus requirement for first natural frequency
figure;
hold on;
plot([0, best_n], [best_h, best_h], '--k', 'LineWidth', 1.3); % horizontal
plot([best_n, best_n], [0, best_h], '--k', 'LineWidth', 1.3); % vertical
contourf(N, H, freqN_H, [threshold, max(freqN_H(:))], 'FaceColor', 'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
greenstar = plot(best_n, best_h, 'p', 'MarkerSize', 12, 'MarkerEdgeColor', 'g', ...
     'MarkerFaceColor', 'g');
grid on;
% Labels and Title
xlabel('n [-]', 'Interpreter', 'latex');
ylabel('h [m]', 'Interpreter', 'latex');
% Create a hidden patch for the legend (same color & transparency)
ylim([0.1, 10])
hPatch = fill(NaN, NaN, 'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add legend
legend([hPatch, greenstar], ...
       {'First natural frequency $>$ 50 Hz', ...
        sprintf('Optimal point (%.2f, %.2f)', best_n, best_h)}, ...
       'Location', 'southeast', 'Interpreter', 'Latex', 'Fontsize', 12);

% Frequency plot
figure;
surf(N, H, freqN_H); hold on
greenstar = plot3(best_n, best_h, freqN_H(row, col)+1, 'p', 'MarkerSize', 12, 'MarkerEdgeColor', 'g', ...
     'MarkerFaceColor', 'g');
xlabel('n [-]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('h [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
zlabel('Frequency [Hz]', 'Interpreter', 'latex');
%title('3D Plot of Frequency vs n and h');
colormap(jet)
h = colorbar; % Adds color bar
xlabel(h, '[Hz]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter','Latex');
h.Label.Rotation = 0; % Set to 0 for horizontal text
shading interp; % Optional: Interpolates shading for smoother surface

% Volume plot
figure;
surf(N, H, VolumeN_H); hold on
greenstar = plot3(best_n, best_h, VolumeN_H(row, col)+0.05, 'p', 'MarkerSize', 12, 'MarkerEdgeColor', 'g', ...
     'MarkerFaceColor', 'g');
xlabel('n [-]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('h [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
zlabel('Volume [m${}^3$]', 'Interpreter', 'latex');
%title('3D Plot of Frequency vs n and h');
colormap(jet)
h = colorbar; % Adds color bar
xlabel(h, '[m${}^3$]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter','Latex');
h.Label.Rotation = 0; % Set to 0 for horizontal text
shading interp; % Optional: Interpolates shading for smoother surface


% Max stress plot
figure;
surf(N, H, maxStressN_H); hold on
greenstar = plot3(best_n, best_h, maxStressN_H(row, col)+5, 'p', 'MarkerSize', 12, 'MarkerEdgeColor', 'g', ...
     'MarkerFaceColor', 'g');
xlabel('n [-]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('h [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
zlabel('Max axial stress [MPa]', 'Interpreter', 'latex');
%title('3D Plot of Frequency vs n and h');
colormap(jet)
h = colorbar; % Adds color bar
xlabel(h, '[MPa]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter','Latex');
%h.Label.Rotation = 0; % Set to 0 for horizontal text
shading interp; % Optional: Interpolates shading for smoother surface

% Min stress plot
figure;
surf(N, H, minStressN_H); hold on
greenstar = plot3(best_n, best_h, minStressN_H(row, col)+10, 'p', 'MarkerSize', 12, 'MarkerEdgeColor', 'g', ...
     'MarkerFaceColor', 'g');
xlabel('n [-]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('h [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
zlabel('Min axial stress [MPa]', 'Interpreter', 'latex');
%title('3D Plot of Frequency vs n and h');
colormap(jet)
h = colorbar; % Adds color bar
xlabel(h, '[MPa]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter','Latex');
%h.Label.Rotation = 0; % Set to 0 for horizontal text
shading interp; % Optional: Interpolates shading for smoother surface