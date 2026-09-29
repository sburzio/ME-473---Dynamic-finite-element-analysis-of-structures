% Mini-Project 3 Code: Natural frequencies and mode shapes of a square plate

% ME-473: Dynamic finite element analysis of structures

% Group F:
% Adrien MAITROT
% Gaël GERBER
% Marko MITRIC
% Hippolyte SOULIER

% Final version done June 1st 2025

clc; clear; close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters

E   = 210e9;        % [Pa]      -   Young's modulus
nu  = 0.3;          % [-]       -   Poisson's ratio
rho = 7850;         % [kg/m^3]  -   Density
h   = 0.005;        % [m]       -   Thickness
L   = 0.5;          % [m]       -   Length

mesh_sizes = [2, 4, 8, 16];     % n x n meshes

% Fréquences analytiques pour modes (m,n) classiques

modes = [1 1; 2 1; 1 2; 2 2];
n_modes = 4;
D = E*h^3/(12*(1 - nu^2));
freq_analytical = zeros(n_modes,1);

for i=1:n_modes
    m_para = modes(i,1);
    n_para = modes(i,2);
    freq_analytical(i) = (pi/2) * sqrt(D/(rho*h)) * (m_para^2 + n_para^2) / (L^2);
end
    

%% 1 & 2: Frequency variation from mesh size using AMC

AMC_All_Modes   = cell(1,length(mesh_sizes));
AMC_All_Freq    = cell(1,length(mesh_sizes));
AMC_All_nodes   = cell(1,length(mesh_sizes));
AMC_All_errors  = zeros(4,length(mesh_sizes));

iter = 1;

fprintf('--- Question 1 & 2 : Frequency for different mesh ---\n');

for n = mesh_sizes

    fprintf('Mesh for %dx%d AMC elements \n', n, n);
    
    % Building rectangular mesh

    [nodes, connectivity] = createRectangularMesh(n,n,L,L);

    % Building stiffness and mass matrices using AMC elements

    [K_AMC, M_AMC] = buildMatricesAMC(nodes, connectivity, E, nu, rho, h);
    
    % Getting the free degrees of freedom 

    AMC_FreeDoFs = getFreeDofs_AMC(n);

    % Computing eigenvalues and eigenvector to deduce the natural
    % frequencies and modeshapes

    [freq_AMC, ModeShape_AMC] = funDynamicAnalysis(K_AMC, M_AMC, AMC_FreeDoFs);

    AMC_All_Modes{iter} = ModeShape_AMC;
    AMC_All_Freq{iter} = freq_AMC;
    AMC_All_nodes{iter} = nodes;
    AMC_All_errors(iter,:) = 100*(freq_AMC - freq_analytical')./(freq_analytical');
    
    % Printing results
    fprintf('Mode | Analytical (Hz) | AMC (Hz)\n');
    for i=1:n_modes
        fprintf('%4d | %15.2f | %8.2f\n', i, freq_analytical(i), freq_AMC(i));
    end
    fprintf('\n');

    iter = iter+1;
end


%% 3: Frequency variation from mesh size using AMC

ratios = 10:10:100;
n_ratios = length(ratios);
freq_thin_h   = zeros(1,n_ratios);
freq_thick_h  = zeros(1,n_ratios);
CR_All_Modes   = cell(1,n_ratios);
CR_All_Freq    = cell(1,n_ratios);
CR_All_nodes   = cell(1,n_ratios);

% Mesh size nxn
n_3 = 8;

fprintf('--- 3 : Varying the thickness ratio ---\n\n');
fprintf('1st mode analysis %dx%d CR elements \n\n', n_3, n_3);
fprintf('Ratio l/h | Thin plate (Hz) | Thick plate (Hz) | CR (Hz)\n');


for i = 1:n_ratios
    
    % Building rectangular mesh

    [nodes, connectivity] = createRectangularMesh(n_3,n_3,L,L);

    % Building stiffness and mass matrices using AMC elements

    h_ratios = L/ratios(i);

    D_h = E*h_ratios^3/(12*(1 - nu^2));

    [K_CR, M_CR] = buildMatricesCR(nodes, connectivity, E, nu, rho, h_ratios);
    
    % Getting the free degrees of freedom 

    CR_FreeDoFs = getFreeDofs_CR(n_3);

    % Computing eigenvalues and eigenvector to deduce the natural
    % frequencies and modeshapes

    [freq_CR, ModeShape_CR] = funDynamicAnalysis(K_CR, M_CR, CR_FreeDoFs);

    CR_All_Modes{i} = ModeShape_CR;
    CR_All_Freq{i} = freq_CR;
    CR_All_nodes{i} = nodes;

    % Computing analytical frequency

    freq_thin_h(i)  = pi/(L^2) * sqrt(D_h/(rho*h_ratios));
    freq_thick_h(i) = computeThickFrequency(L, L, h_ratios, nu, E, rho, 5/6, 1, 1);
    
    % Printing results
    fprintf('%9d | %15.2f | %16.2f | %8.2f\n', ratios(i), freq_thin_h(i), freq_thick_h(i), CR_All_Freq{i}(1));
end

fprintf('\n');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots section

% Plots for questions 1 and 2, mode shapes for different mesh sizes

for mesh = 1:length(mesh_sizes)
    plotModeShape2D(AMC_All_nodes{mesh}, AMC_All_Modes{mesh}(:,1), 3, 1, 1)
    plotModeShape2D(AMC_All_nodes{mesh}, AMC_All_Modes{mesh}(:,2), 3, 1, 2)
    plotModeShape2D(AMC_All_nodes{mesh}, AMC_All_Modes{mesh}(:,3), 3, 2, 1)
    plotModeShape2D(AMC_All_nodes{mesh}, AMC_All_Modes{mesh}(:,4), 3, 2, 2)
end


% Small animation of the 16x16 mesh for mode(2,2)

animateModeShape3D(AMC_All_nodes{4}, AMC_All_Modes{4}(:,4), 3)


% Plot of the error between the analytical and the numerical solution by
% varying mesh size

figure()
hold on
plot([2, 4, 8, 16], AMC_All_errors(:,1), 'o-', 'LineWidth', 1.5)
plot([2, 4, 8, 16], AMC_All_errors(:,2), 'o-', 'LineWidth', 1.5)
plot([2, 4, 8, 16], AMC_All_errors(:,3), 'o-', 'LineWidth', 1.5)
plot([2, 4, 8, 16], AMC_All_errors(:,4), 'o-', 'LineWidth', 1.5)
xlabel('Number of elements [-]')
ylabel('Error [%]')
legend('Mode(1,1)', 'Mode(1,2)', 'Mode(2,1)', 'Mode(2,2)')
grid on


% Plot of the comparison between thick and thin plate theory

CR_f11 = zeros(n_ratios,1);

for i = 1:n_ratios
    CR_f11(i) = CR_All_Freq{i}(1);
end

figure()
hold on
plot(ratios, CR_f11, 'o-', 'LineWidth', 1.5)
plot(ratios, freq_thick_h, 'o-', 'LineWidth', 1.5)
plot(ratios, freq_thin_h, 'o-', 'LineWidth', 1.5)
xlabel('Geometrical ratio l/h [-]')
ylabel('Frequency [Hz]')
legend('FEA using CR', 'Thick plate', 'Exact')
grid on


% Plot of l/h = 10 and l/h = 100 second modeshape

plotModeShape2D(CR_All_nodes{1}, CR_All_Modes{1}(:,2), 4, 1, 2)
plotModeShape2D(CR_All_nodes{10}, CR_All_Modes{10}(:,2), 4, 1, 2)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

function [nodes, connectivity] = createRectangularMesh(n_elem_x,n_elem_y, a, b)
    % Create a rectangular mesh with 
    % number_of_elements_X * number_of_elements_Y  
    % elements, over a rectangle of size a x b.
    %
    % Output:
    % nodes: (N_nodes x 2) array with (x,y) coordinates
    % connectivity: (N_elements x 4) array with node indices (quadrilateral elements)

    % Number of nodes in each direction
    n_nodes_x = n_elem_x + 1;
    n_nodes_y = n_elem_y + 1;

    % Generate grid points
    x = linspace(0, a, n_nodes_x);
    y = linspace(0, b, n_nodes_y);
    [X, Y] = meshgrid(x, y);

    % Node list: flatten the (X, Y) matrices
    nodes = [X(:), Y(:)];

    % Connectivity: define each element by 4 nodes
    connectivity = cell(1, n_elem_x * n_elem_y); 

    % connectivity = zeros(n_elem_x * n_elem_y, 4);
    elem = 1;
    for i = 1:n_elem_x
        for j = 1:n_elem_y
            n1 = (i-1)*n_nodes_y + j;       % bottom-left
            n2 = n1 + n_nodes_y;            % bottom-right
            n3 = n2 + 1;                    % top-right
            n4 = n1 + 1;                    % top-left
            connectivity{elem} = [n1, n2, n3, n4];
            elem = elem + 1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [K, M] = buildMatricesAMC(nodes, connectivity, E, nu, rho, h)

    number_elements = size(connectivity,2);
    number_Nodes = size(nodes,1);
    number_dofs = number_Nodes*3; 

    eL = sqrt((nodes(2,1)-nodes(1,1))^2 + (nodes(2,2)-nodes(1,2))^2);

    C = E/(1-nu^2)*[1,nu,0; nu,1,0; 0,0,(1-nu)/2];

    % Computing the element stiffness and mass matrix

    xi1 = sym('xi1');
    xi2 = sym('xi2');

    ah_i = @(xi1_i, xi2_i) [(1 + xi1_i*xi1)*(1 + xi2_i*xi2)*(2 + xi1_i*xi1 + xi2_i*xi2 - xi1^2 - xi2^2),...
                            eL/16 * (1 + xi1_i*xi1)*(xi2_i + xi2)*(xi2^2-1),...
                            -eL/16 * (1 + xi2_i*xi2)*(xi1_i + xi1)*(xi1^2-1)];

    aH = [ah_i(-1,-1), ah_i(1,-1), ah_i(1,1), ah_i(-1,1)];

    aB = [4/eL^2*diff(diff(aH,xi1),xi1); 4/eL^2*diff(diff(aH,xi2),xi2); 8/eL^2*diff(diff(aH,xi2),xi1)];

    eK = eL^2*h^3/12*int(int(aB'*C*aB,xi2,([-1,1])),xi1,[-1,1]);

    eM = eL^2*int(int(aH'*rho*h*aH,xi2,([-1,1])),xi1,[-1,1]);
    
    % Assembling the matrices

    K = zeros(number_dofs);
    M = zeros(number_dofs);

    for e=1:number_elements
        conn_elem = connectivity{e};  
        indices_dofs = reshape([3*conn_elem-2; 3*conn_elem-1; 3*conn_elem], [], 1);
    
        K(indices_dofs,indices_dofs) = K(indices_dofs,indices_dofs) + eK;
        M(indices_dofs,indices_dofs) = M(indices_dofs,indices_dofs) + eM;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [K, M] = buildMatricesCR(nodes, connectivity, E, nu, rho, h)

    number_elements = size(connectivity,2);
    number_Nodes = size(nodes,1);
    number_dofs = number_Nodes*4; 

    eL = sqrt((nodes(2,1)-nodes(1,1))^2 + (nodes(2,2)-nodes(1,2))^2);

    C = E/(1-nu^2)*[1,nu,0; nu,1,0; 0,0,(1-nu)/2];

    % Computing the element stiffness and mass matrix

    xi1 = sym('xi1');
    xi2 = sym('xi2');

    fi = @(xi_i, xi) (-xi_i*xi^3 + 3*xi_i*xi + 2)/4;
    gi= @(xi_i, xi)  (xi^3 + xi_i*xi^2 - xi - xi_i)/4;

    ah_i = @(xi1_i, xi2_i)[fi(xi1_i,xi1)*fi(xi2_i,xi2), ...
                            eL/2 * fi(xi1_i,xi1)*gi(xi2_i,xi2), ...
                            -eL/2 * gi(xi1_i,xi1)*fi(xi2_i,xi2), ...
                            eL^2/4 * gi(xi1_i,xi1)*gi(xi2_i,xi2)];

    aH = [ah_i(-1,-1), ah_i(1,-1), ah_i(1,1), ah_i(-1,1)];

    aB = [4/eL^2*diff(diff(aH,xi1),xi1); 4/eL^2*diff(diff(aH,xi2),xi2); 8/eL^2*diff(diff(aH,xi2),xi1)];

    eK = eL^2*h^3/12*int(int(aB'*C*aB,xi2,([-1,1])),xi1,[-1,1]);

    eM = eL^2*int(int(aH'*rho*h*aH,xi2,([-1,1])),xi1,[-1,1]);
    
    % Assembling the matrices

    K = zeros(number_dofs);
    M = zeros(number_dofs);

    for e=1:number_elements
        conn_elem = connectivity{e};  
        indices_dofs = reshape([4*conn_elem-3; 4*conn_elem-2; 4*conn_elem-1; 4*conn_elem], [], 1);
    
        K(indices_dofs,indices_dofs) = K(indices_dofs,indices_dofs) + eK;
        M(indices_dofs,indices_dofs) = M(indices_dofs,indices_dofs) + eM;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [free_dofs] = getFreeDofs_AMC(n)

    % Finding corner nodes

    bottom_left  =  1;
    top_left     =  n + 1;
    bottom_right = (n + 1) * n + 1;
    top_right    = (n + 1) * (n + 1);

    corner_nodes_indices = [bottom_left, bottom_right, top_right, top_left];

    % Finding edge nodes

    left_edge   = (bottom_left + 1):(top_left - 1);
    right_edge  = (bottom_right + 1):(top_right - 1);
    bottom_edge = (top_left + 1):(n + 1):(bottom_right - 1);
    top_edge    = (top_left + n + 1):(n + 1):(top_right - 1);

    edge_X_indices = [bottom_edge, top_edge]; % edges parallel to X
    edge_Y_indices = [right_edge, left_edge]; % edges parallel to Y

    % Defining the constrained DoFs

    constrained_dofs = []; 

    for i = 1:length(corner_nodes_indices)
        node_index = corner_nodes_indices(i);    
        constrained_dofs = [constrained_dofs, 3*node_index-2, 3*node_index-1, 3*node_index];
    end

    for i = 1:length(edge_X_indices)
        node_index = edge_X_indices(i);       
        constrained_dofs = [constrained_dofs, 3*node_index-2, 3*node_index];
    end

    for i = 1:length(edge_Y_indices)
        node_index = edge_Y_indices(i);       
        constrained_dofs = [constrained_dofs, 3*node_index-2, 3*node_index-1];
    end

    constrained_dofs = sort(constrained_dofs);
    
    free_dofs = setdiff(1:(3*top_right), constrained_dofs);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [free_dofs] = getFreeDofs_CR(n)

    % Finding corner nodes

    bottom_left  =  1;
    top_left     =  n + 1;
    bottom_right = (n + 1) * n + 1;
    top_right    = (n + 1) * (n + 1);

    corner_nodes_indices = [bottom_left, bottom_right, top_right, top_left];

    % Finding edge nodes

    left_edge   = (bottom_left + 1):(top_left - 1);
    right_edge  = (bottom_right + 1):(top_right - 1);
    bottom_edge = (top_left + 1):(n + 1):(bottom_right - 1);
    top_edge    = (top_left + n + 1):(n + 1):(top_right - 1);

    edge_X_indices = [bottom_edge, top_edge]; % edges parallel to X
    edge_Y_indices = [right_edge, left_edge]; % edges parallel to Y

    % Defining the constrained DoFs

    constrained_dofs = []; 

    for i = 1:length(corner_nodes_indices)
        node_index = corner_nodes_indices(i);    
        constrained_dofs = [constrained_dofs, 4*node_index-3, 4*node_index-2, 4*node_index-1, 4*node_index];
    end

    for i = 1:length(edge_X_indices)
        node_index = edge_X_indices(i);       
        constrained_dofs = [constrained_dofs, 4*node_index-3, 4*node_index-1, 4*node_index];
    end

    for i = 1:length(edge_Y_indices)
        node_index = edge_Y_indices(i);       
        constrained_dofs = [constrained_dofs, 4*node_index-3, 4*node_index-2, 4*node_index];
    end

    constrained_dofs = sort(constrained_dofs);
    
    free_dofs = setdiff(1:(4*top_right), constrained_dofs);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% The function below computes the dynamical analysis (frequency and mode shapes)

function [Freq, ModeShape] = funDynamicAnalysis(Stiffness, Mass, ActiveDoF)

    % Computing the eigen vectors and values

    [EigVectors, EigValues] = eig(Mass(ActiveDoF, ActiveDoF)\Stiffness(ActiveDoF, ActiveDoF));

    EigValues = sqrt(EigValues);

    % Finding the 4 smallest natural frequencies in Hertz 

    Freq = mink(maxk(EigValues,1,1),4)./(2*pi); 

    % Computing the Mode Shapes

    ModeShape = zeros(size(Mass,1),4);

    for j = 1:4

         % Finding the location of the frequency in eigen values matrix

        [~, col] = find(EigValues./(2*pi) == Freq(j));

        % Replacing the active DoF modal shape

        ModeShape(ActiveDoF,j) = EigVectors(:,col);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function animateModeShape3D(nodes, ModeShape, nodeDOF)

    % Extract vertical displacement (assuming w is the first DoF)
    modeIndex = 1:nodeDOF:(nodeDOF*size(nodes,1));
    w = ModeShape(modeIndex);

    % Normalize displacement
    scaleFactor = 0.1 / max(abs(w));
    w_scaled = w * scaleFactor;

    % Scattered interpolant
    F = scatteredInterpolant(nodes(:,1), nodes(:,2), w_scaled, 'natural', 'none');

    % Grid for interpolation
    nGrid = 100;
    xq = linspace(min(nodes(:,1)), max(nodes(:,1)), nGrid);
    yq = linspace(min(nodes(:,2)), max(nodes(:,2)), nGrid);
    [Xq, Yq] = meshgrid(xq, yq);
    Zq_base = F(Xq, Yq);

    % Animation parameters
    nFrames = 200;
    omega = pi/2;  % arbitrary frequency
    T = linspace(0, 2*pi, nFrames);  % simulate one full period
    maxAmp = max(abs(Zq_base), [], 'all');

    % Setup figure
    figure;
    hSurf = surf(Xq, Yq, F(Xq, Yq), 'EdgeColor', 'none');
    colormap jet;
    colorbar;
    clim([-maxAmp, maxAmp]); % Fix colormap
    zlim([-maxAmp, maxAmp]);  % Fix Z axis
    xlim([min(xq), max(xq)]);
    ylim([min(yq), max(yq)]);
    axis equal;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Amplitude [-]');
    view(3);

    % Animation loop
    for t = T
        Zq = F(Xq, Yq) * sin(omega * t);
        set(hSurf, 'ZData', Zq);
        title(sprintf('Vibrating Mode Shape'));
        drawnow;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function plotModeShape2D(nodes, ModeShape, nodeDOF, modex, modey)

    % Getting the DoFs for the displacement (assuming vertical is first)
    modeIndex = 1:nodeDOF:(nodeDOF*size(nodes,1));

    % Extract vertical displacement vector 
    w = ModeShape(modeIndex);

    % Normalize displacement
    max_index = (abs(w) == max(abs(w)));
    scaleFactor = 0.1 / w(max_index);
    w_scaled = w * scaleFactor;

    % Create scattered interpolant for smooth coloring
    F = scatteredInterpolant(nodes(:,1), nodes(:,2), w_scaled, 'natural', 'none');

    % Create a regular grid for visualization
    nGrid = 300;
    xq = linspace(min(nodes(:,1)), max(nodes(:,1)), nGrid);
    yq = linspace(min(nodes(:,2)), max(nodes(:,2)), nGrid);
    [Xq, Yq] = meshgrid(xq, yq);
    Zq = F(Xq, Yq);

    % Plot as a 2D colormap
    figure;
    contourf(Xq, Yq, Zq, 30, 'LineColor', 'none');
    axis equal;
    axis tight;
    colormap jet;
    colorbar;
    xlabel('X [m]');
    ylabel('Y [m]');
    titlegraph = "2D Mode Shape (Amplitude colormap) - mode(" + num2str(modex) + ", " + num2str(modey) + ")";
    title(titlegraph);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function omega_thick = computeThickFrequency(length_X, length_Y, h, nu, E, rho, k, n, m)
    % Computes the exact natural frequency of indices n m 
    % for a plate using given parameters.
    % Shear modulus
    G = E / (2 * (1 + nu));
    % Flexural rigidity
    D = E * h^3 / (12 * (1 - nu^2));
    % Wavenumbers
    alpha = m * pi / length_X;
    beta = n * pi / length_Y;
    
    % Matrix of flexural rigidities
    D11 = D;
    D12 = D * nu;
    D22 = D;
    D66 = D * (1 - nu)/2;
    % Matrix of shear rigidities
    S44 = h * k * G;
    S55 = h * k * G;
    % Terms computation
    P11 = D11 * alpha^2 + D66 * beta^2 + S55;
    P12 = (D12 + D66) * alpha * beta;
    P13 = S55 * alpha;
    P22 = D66 * alpha^2 + D22 * beta^2 + S44;
    P23 = S44 * beta;
    P33 = S55 * alpha^2 + S44 * beta^2;
    % Intermediate determinant
    Q = P11 * P22 - P12^2;
    % Frequency squared
    omega_squared = (Q * P33 + 2 * P12 * P23 * P13 - P22 * P13^2 - P11 * P23^2) / (rho * h * Q);
    % Natural frequency
    omega_thick = sqrt(omega_squared)/(2*pi);
end