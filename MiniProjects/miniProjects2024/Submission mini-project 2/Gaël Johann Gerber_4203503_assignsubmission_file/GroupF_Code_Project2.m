% Mini-Project 2 Code: Design of a 2d truss bridge

% ME-473: Dynamic finite element analysis of structures

% Group F:
% Adrien MAITROT
% Gaël GERBER
% Marko MITRIC
% Hippolyte SOULIER

% Final version done 15th April 2025

clc; clear; close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Parameters (non-variable)

E = [200e9, 69e9, 12e9];    % [Pa]      -   Young's modulus of steel, alu and wood
rho = [7850, 2700, 600];    % [kg/m^3]  -   Steel Density
p = 50000;                  % [N]       -   Applied external load

anchorCoordinates = [0 0;...
                    -1 0;...
                    20 0;...
                    21 0];  % [m]       -   Anchors x and y coordinates (nodes 1,2,3 and 4)

PrescribedDoF = [1; 2;...
                 3; 4;...
                 5; 6;...
                 7; 8];     % No displacement on anchor degrees of freedom

% Design requirements:

u_max = 0.04;               % [m]       -   Maximal node displacement
V_max = [1, 1.7, 2];        % [m^3]     -   Volume limit for steel, alu and wood
f_min = 1.5;                % [Hz]      -   Frequency limit
f_max = 10;                 % [Hz]      -   Frequency limit
yield = [350, 110, 30]*10^6;% [MPa]     -   Yield stress for steel, alu and wood


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plots
modeplot = true;





%% Steel Bridge

%%%%%%%%%%%%%%%%%%%%%%%    Geometric Construction    %%%%%%%%%%%%%%%%%%%%%%

% Defining the Nodes location and connectivity table

theta = linspace(120,60,15);

SteelNodes = [anchorCoordinates; (10.5*(cosd(theta')/sind(30)+1)-0.5) (10.5*(sind(theta')/sind(30)-1/tand(30))+0.5); 10 2];
SteelConnec = [1 2; 3 4; 1 5; 2 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11; 11 12;...
                12 13; 13 14; 14 15; 15 16; 16 17; 17 18; 18 19; 19 3; 19 4;...
                1 20; 3 20; 5 20 ; 6 20; 7 20; 8 20; 9 20; 10 20;11 20; ...
                12 20; 13 20; 14 20; 15 20; 16 20; 17 20; 18 20; 19 20];

SteelDoF = 2*size(SteelNodes,1);

% Defining the Area of each trusses, if we want to have different areas

A_Steel = 15*10^(-4)*ones(size(SteelNodes,1));  % [m^2]

% Plotting the bridge 

DrawBridge(SteelNodes, SteelConnec, 'Title', 'Steel Bridge Structure')

% Computing the stiffness and mass matrix and the total volume

[K_Steel, M_Steel, VTot_Steel] = funMatrices(SteelDoF, SteelConnec, SteelNodes, E(1), A_Steel, rho(1));

% Applying Boundary Conditions

SteelActiveDof = setdiff(transpose((1:SteelDoF)), PrescribedDoF);

SteelExtLoad = zeros(SteelDoF,1);
SteelExtLoad(12*2) = -p;


%%%%%%%%%%%%%%%%%%%%%%%%%%    Static Analysis    %%%%%%%%%%%%%%%%%%%%%%%%%%

[U_Steel, U_max_Steel, Node_max_Steel, stress_Steel, strain_Steel] = funStaticAnalysis(SteelExtLoad, K_Steel, SteelActiveDof, size(SteelNodes,1), SteelNodes, SteelConnec, E(1));

SteelNodesDef = SteelNodes + [U_Steel(1:2:end), U_Steel(2:2:end)];

DrawBridge(SteelNodesDef, SteelConnec, 'Title', 'Steel Bridge Deformation and stress', 'ElementStress', stress_Steel, 'ForceNode', 12, 'ShowElementNumbers', false)


%%%%%%%%%%%%%%%%%%%%%%%%%%    Dynamic Analysis   %%%%%%%%%%%%%%%%%%%%%%%%%%

[Freq_Steel, ModeShape_Steel] = funDynamicAnalysis(K_Steel, M_Steel, SteelActiveDof);

SteelNodesMode = zeros(size(SteelNodes,1),2,4);

% Defining a scale factor for the modeshape

Steelscalefactor = [3,3,3,3];

for j = 1:4
    SteelNodesMode(:,:,j) = SteelNodes + Steelscalefactor(j)*[ModeShape_Steel(1:2:end,j), ModeShape_Steel(2:2:end,j)];

    % Plotting the modeshape

    if modeplot
        DrawBridge(SteelNodesMode(:,:,j), SteelConnec, 'Title', "Steel Bridge Modeshape " + int2str(j))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%        Results        %%%%%%%%%%%%%%%%%%%%%%%%%%

% Display Criteria Results

fprintf('------- Steel Bridge -------\n')

if(max(abs(stress_Steel)) > yield(1))
    fprintf('\nExceeded stress limit: the yield stress is %.0f MPa, actual design has %.2f MPa', yield(1)*10^(-6),max(abs(stress_Steel))*10^(-6))
else
    fprintf('\nThe maximal stress is smaller than %.0f MPa, the actual design has %.2f MPa,', yield(1)*10^(-6),max(abs(stress_Steel))*10^(-6))
end

if(U_max_Steel > u_max)
    fprintf('\nExceeded displacement limit: the nodal displacement must be smaller than 40 mm, actual design has %.3f mm', U_max_Steel*1000)
else
    fprintf('\nThe maximal nodal displacement is smaller than 40 mm, the actual design has %.3f mm', U_max_Steel*1000)
end

if(VTot_Steel > V_max(1))
    fprintf('\nExceeded volume limit: the total volume must be smaller than 1 m^3, actual design is %.3f m^3', VTot_Steel)
else
    fprintf('\nThe total volume is smaller than 1 m^3, actual design is %.3f m^3', VTot_Steel)
end

if((Freq_Steel(1) >= f_min) && (Freq_Steel(1) <= f_max))
    fprintf('\nWarning: the first natural frequency is within range, actual design has %.3f Hz', Freq_Steel(1))
else
    fprintf('\nThe first natural frequency is outside range, actual design has %.3f Hz', Freq_Steel(1))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Aluminum Bridge

%%%%%%%%%%%%%%%%%%%%%%%    Geometric Construction    %%%%%%%%%%%%%%%%%%%%%%

% Defining the Nodes location and connectivity table
AluNodes = [anchorCoordinates; 10 4; -1 2; 0 2; -1 4; 0 4;
     -1 6; 0 6;  -1 8; 0 8;  20 2; 21 2;  20 4; 21 4;
     20 6; 21 6;  20 8; 21 8; 2.5 6.5; 5 5.7; 7.5 5.3;
    10 5; 12.5 5.3; 15 5.7; 17.5 6.5; 2.5 4; 5 4; 7.5 4; 12.5 4; 15 4; 17.5 4;];
AluConnec = [1 2; 3 4;
    1 6; 2 7; 6 7; 7 8; 6 9; 8 9; 8 11; 9 10; 10 11; 10 13; 11 12; 12 13;
    2 6; 6 8; 8 10; 10 12; 1 7; 7 9; 9 11; 11 13;
    3 15; 4 14; 14 15; 14 17; 15 16; 16 17; 16 19; 17 18; 18 19; 18 21; 19 20; 20 21;
    3 14; 14 16; 16 18; 18 20; 4 15; 15 17; 17 19; 19 21;
    9 29; 29 30; 30 31; 31 5; 5 32; 32 33; 33 34; 34 16;
    13 22; 22 23; 23 24; 24 25; 25 26; 26 27; 27 28; 28 20;
    13 29; 9 22; 22 30; 29 23; 23 31; 30 24; 24 5; 31 25; 25 32; 
    5 26; 26 33; 32 27; 33 28; 28 16; 27 34; 34 20;
    1 5; 2 31; 7 30; 6 29; 3 5; 4 32; 14 33; 15 34; 1 30; 1 31; 3 32; 3 33; 2 30; 2 5; 4 33; 4 5];

AluDoF = 2*size(AluNodes,1);

% Defining the Area of each trusses

A_Alu = 20*10^(-4)*ones(size(AluNodes,1));  % [m^2]

% Plotting the bridge 

DrawBridge(AluNodes, AluConnec, 'Title', 'Aluminum Bridge Structure')

% Computing the stiffness and mass matrix and the total volume

[K_Alu, M_Alu, VTot_Alu] = funMatrices(AluDoF, AluConnec, AluNodes, E(2), A_Alu, rho(2));

% Applying Boundary Conditions

AluActiveDof = setdiff(transpose((1:AluDoF)), PrescribedDoF);

AluExtLoad = zeros(AluDoF,1);
AluExtLoad(5*2) = -p;


%%%%%%%%%%%%%%%%%%%%%%%%%%    Static Analysis    %%%%%%%%%%%%%%%%%%%%%%%%%%

[U_Alu, U_max_Alu, Node_max_Alu, stress_Alu, strain_Alu] = funStaticAnalysis(AluExtLoad, K_Alu, AluActiveDof, size(AluNodes,1), AluNodes, AluConnec, E(2));

AluNodesDef = AluNodes + [U_Alu(1:2:end), U_Alu(2:2:end)];

DrawBridge(AluNodesDef, AluConnec, 'Title', 'Aluminum Bridge Deformation and stress', 'ElementStress', stress_Alu, 'ForceNode', 5, 'ShowElementNumbers', false)


%%%%%%%%%%%%%%%%%%%%%%%%%%    Dynamic Analysis   %%%%%%%%%%%%%%%%%%%%%%%%%%

[Freq_Alu, ModeShape_Alu] = funDynamicAnalysis(K_Alu, M_Alu, AluActiveDof);

AluNodesMode = zeros(size(AluNodes,1),2,4);

% Defining a scale factor for the modeshape

Aluscalefactor = [3,3,3,3];

for j = 1:4
    AluNodesMode(:,:,j) = AluNodes + Aluscalefactor(j)*[ModeShape_Alu(1:2:end,j), ModeShape_Alu(2:2:end,j)];

    % Plotting the modeshape

    if modeplot
        DrawBridge(AluNodesMode(:,:,j), AluConnec, 'Title', "Aluminum Bridge Modeshape " + int2str(j))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%        Results        %%%%%%%%%%%%%%%%%%%%%%%%%%

% Display Criteria Results

fprintf('\n\n------- Alu Bridge -------\n')

if(max(abs(stress_Alu)) > yield(2))
    fprintf('\nExceeded stress limit: the yield stress is %.0f MPa, actual design has %.2f MPa', yield(2)*10^(-6),max(abs(stress_Alu))*10^(-6))
else
    fprintf('\nThe maximal stress is smaller than %.0f MPa, the actual design has %.2f MPa,', yield(2)*10^(-6),max(abs(stress_Alu))*10^(-6))
end

if(U_max_Alu > u_max)
    fprintf('\nExceeded displacement limit: the nodal displacement must be smaller than 40 mm, actual design has %.3f mm', U_max_Alu*1000)
else
    fprintf('\nThe maximal nodal displacement is smaller than 40 mm, the actual design has %.3f mm', U_max_Alu*1000)
end

if(VTot_Alu > V_max(2))
    fprintf('\nExceeded volume limit: the total volume must be smaller than 1.7 m^3, actual design is %.3f m^3', VTot_Alu)
else
    fprintf('\nThe total volume is smaller than 1.7 m^3, actual design is %.3f m^3', VTot_Alu)
end


if((Freq_Alu(1) >= f_min) && (Freq_Alu(1) <= f_max))
    fprintf('\nWarning: the first natural frequency is within range, actual design has %.3f Hz', Freq_Alu(1))
else
    fprintf('\nThe first natural frequency is outside range, actual design has %.3f Hz ', Freq_Alu(1))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Wood Bridge

%%%%%%%%%%%%%%%%%%%%%%%    Geometric Construction    %%%%%%%%%%%%%%%%%%%%%%

% Defining the Nodes location and connectivity table

WoodNodes = [anchorCoordinates; 10 7; 10 2; 2 4; 18 4];
WoodConnec = [1 2; 3 4; 1 7; 2 7; 3 8; 4 8; 1 6; 3 6; 5 6; 5 7; 5 8; 1 5; 7 6; 6 8; 3 5];

WoodDoF = 2*size(WoodNodes,1);

% Defining the Area of each trusses

A_Wood = 30*10^(-4)*ones(size(WoodNodes,1));  % [m^2]

% Plotting the bridge 

DrawBridge(WoodNodes, WoodConnec, 'Title', 'Wood Bridge Structure')

% Computing the stiffness and mass matrix and the total volume

[K_Wood, M_Wood, VTot_Wood] = funMatrices(WoodDoF, WoodConnec, WoodNodes, E(3), A_Wood, rho(3));

% Applying Boundary Conditions

WoodActiveDof = setdiff(transpose((1:WoodDoF)), PrescribedDoF);

WoodExtLoad = zeros(WoodDoF,1);
WoodExtLoad(5*2) = -p;


%%%%%%%%%%%%%%%%%%%%%%%%%%    Static Analysis    %%%%%%%%%%%%%%%%%%%%%%%%%%

[U_Wood, U_max_Wood, Node_max_Wood, stress_Wood, strain_Wood] = funStaticAnalysis(WoodExtLoad, K_Wood, WoodActiveDof, size(WoodNodes,1), WoodNodes, WoodConnec, E(3));

WoodNodesDef = WoodNodes + [U_Wood(1:2:end), U_Wood(2:2:end)];

DrawBridge(WoodNodesDef, WoodConnec, 'Title', 'Wood Bridge Deformation and stress', 'ElementStress', stress_Wood, 'ForceNode', 5, 'ShowElementNumbers', false)


%%%%%%%%%%%%%%%%%%%%%%%%%%    Dynamic Analysis   %%%%%%%%%%%%%%%%%%%%%%%%%%

[Freq_Wood, ModeShape_Wood] = funDynamicAnalysis(K_Wood, M_Wood, WoodActiveDof);

WoodNodesMode = zeros(size(WoodNodes,1),2,4);

% Defining a scale factor for the modeshape

Woodscalefactor = [3,3,3,3];

for j = 1:4
    WoodNodesMode(:,:,j) = WoodNodes + Woodscalefactor(j)*[ModeShape_Wood(1:2:end,j), ModeShape_Wood(2:2:end,j)];

    % Plotting the modeshape

    if modeplot
        DrawBridge(WoodNodesMode(:,:,j), WoodConnec, 'Title', "Wood Bridge Modeshape " + int2str(j))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%        Results        %%%%%%%%%%%%%%%%%%%%%%%%%%

% Display Criteria Results

fprintf('\n\n------- Wood Bridge -------\n')

if(max(abs(stress_Wood)) > yield(3))
    fprintf('\nExceeded stress limit: the yield stress is %.0f MPa, actual design has %.2f MPa', yield(3)*10^(-6),max(abs(stress_Wood))*10^(-6))
else
    fprintf('\nThe maximal stress is smaller than %.0f MPa, the actual design has %.2f MPa,', yield(3)*10^(-6),max(abs(stress_Wood))*10^(-6))
end

if(U_max_Wood > u_max)
    fprintf('\nExceeded displacement limit: the nodal displacement must be smaller than 40 mm, actual design has %.3f mm', U_max_Wood*1000)
else
    fprintf('\nThe maximal nodal displacement is smaller than 40 mm, the actual design has %.3f mm', U_max_Wood*1000)
end

if(VTot_Wood > V_max(3))
    fprintf('\nExceeded volume limit: the total volume must be smaller than 2 m^3, actual design is %.3f m^3', VTot_Wood)
else
    fprintf('\nThe total volume is smaller than 2 m^3, actual design is %.3f m^3', VTot_Wood)
end

if((Freq_Wood(1) >= f_min) && (Freq_Wood(1) <= f_max))
    fprintf('\nWarning: the first natural frequency is within range, actual design has %.3f Hz\n', Freq_Wood(1))
else
    fprintf('\nThe first natural frequency is outside range, actual design has %.3f Hz\n', Freq_Wood(1))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Functions

% Computation of the assembled stiffness and mass matrices and the total volume

function [Stiffness, TotalMass, VTotal] = funMatrices(GDof, connectivity, nodesCoordinates, E, A, rho)

    % Initializing the assembled matrix and the number of elements

    Stiffness = zeros(GDof);
    TotalMass = zeros(GDof);
    VTotal = 0;
    N_Elements = size(connectivity,1);
    
    % Extracting the X and Y coordinates of the nodes

    Nodes_X = nodesCoordinates(:,1);
    Nodes_Y = nodesCoordinates(:,2);

    for i = 1:N_Elements

        % Extracting the i-th element parameters

        localIndices = connectivity(i,:);
        elementDof = [2*localIndices(1)-1 2*localIndices(1) 2*localIndices(2)-1 2*localIndices(2)] ;
        ElementArea = A(i);
    
        % Computing the x-length, y-length and total length

        xa = Nodes_X(localIndices(2)) - Nodes_X(localIndices(1));
        ya = Nodes_Y(localIndices(2)) - Nodes_Y(localIndices(1));
        L_element = sqrt(xa*xa+ya*ya);

        % Computing the cosine and sine

        C = xa/L_element;
        S = ya/L_element;
    
        % Computing the element stiffness, mass and volume

        ElementStiffness = E*ElementArea/L_element*[C*C C*S -C*C -C*S; C*S S*S -C*S -S*S;-C*C -C*S C*C C*S;-C*S -S*S C*S S*S];
        ElementMass = rho*ElementArea*L_element/6*[2*C*C 2*C*S C*C C*S; 2*C*S 2*S*S C*S S*S; C*C C*S 2*C*C 2*C*S;C*S S*S 2*C*S 2*S*S];
        ElementV = L_element*ElementArea;

        % Assembling the total stifness and mass matrices

        Stiffness(elementDof,elementDof) = Stiffness(elementDof,elementDof) + ElementStiffness;
        TotalMass(elementDof,elementDof) = TotalMass(elementDof,elementDof) + ElementMass;
        VTotal = VTotal + ElementV;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The function below plots the bridge

function DrawBridge(nodesCoordinates, connectivity, varargin)

    % Parse input options
    p = inputParser;
    addParameter(p, 'LineColor', 'b');
    addParameter(p, 'LineWidth', 1);
    addParameter(p, 'MarkerSize', 15);
    addParameter(p, 'ShowNodeNumbers', true);
    addParameter(p, 'Title', '2D Truss Structure');
    addParameter(p, 'FontSize', 8)
    addParameter(p, 'ShowElementNumbers', true)
    addParameter(p, 'ElementStress', []); 
    addParameter(p, 'ForceNode', []);       
    parse(p, varargin{:});
    
    % Extract parsed options
    lineColor = p.Results.LineColor;
    lineWidth = p.Results.LineWidth;
    markerSize = p.Results.MarkerSize;
    showNodeNumbers = p.Results.ShowNodeNumbers;
    fontSize = p.Results.FontSize;
    titleText = p.Results.Title;
    showElementNumbers = p.Results.ShowElementNumbers;
    elementStress = p.Results.ElementStress;
    forceNode = p.Results.ForceNode;

    figure; hold on; axis equal;
    xlabel('X-Coordinate [m]'); ylabel('Y-Coordinate [m]');
    title(titleText, 'Interpreter', 'none');

    % Define colormap 
    cmap = jet(256); % 256-color gradient from blue to red

    if ~isempty(elementStress)
       % Normalize stress values to 1–256
       stressMin = min(elementStress);
       stressMax = max(elementStress);
       stressNorm = round(1 + 255 * (elementStress - stressMin) / (stressMax - stressMin));
    end

    % Plot truss elements
    for i = 1:size(connectivity, 1)
        n1 = connectivity(i, 1);
        n2 = connectivity(i, 2);
        x = [nodesCoordinates(n1, 1), nodesCoordinates(n2, 1)];
        y = [nodesCoordinates(n1, 2), nodesCoordinates(n2, 2)];

        if ~isempty(elementStress)
            colorIdx = stressNorm(i);
            thisColor = cmap(colorIdx, :); % RGB triplet from colormap
        else
            thisColor = lineColor; % fallback
        end

        plot(x, y, '-', 'Color', thisColor, 'LineWidth', lineWidth);

        % Midpoint of element
        midX = mean(x);
        midY = mean(y);

        if showElementNumbers
            % Plot white square with black edge at midpoint
            plot(midX, midY, 'gs', ...
                'MarkerSize', markerSize * 0.8, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'w');

            % Label the truss element number (larger, black font)
            text(midX, midY, sprintf('%d', i), ...
                'Color', 'b', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end

    % Adding the colormap
    if ~isempty(elementStress)
        colormap(jet(256));
        colorbar;
        clim([stressMin, stressMax]);
        ylabel(colorbar, 'Element Stress (Pa)');
    end

    % Plot the anchors 
    for j = 1:4
        
        node_x = nodesCoordinates(j,1);
        node_y = nodesCoordinates(j,2);

        triangle_x = [node_x, node_x-1/2.5, node_x+1/2.5];
        triangle_y = [node_y, node_y-1/1.5, node_y-1/1.5];

        fill(triangle_x, triangle_y, 'w', 'EdgeColor','k', 'LineWidth', 1)
    end


    % Plot nodes with white interior circles
    plot(nodesCoordinates(:,1), nodesCoordinates(:,2), 'o', ...
         'MarkerSize', markerSize*0.5, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

    % Node numbers inside the circle
    if showNodeNumbers
        for i = 1:size(nodesCoordinates, 1)
            text(nodesCoordinates(i,1), nodesCoordinates(i,2), ...
                sprintf('%d', i), ...
                'FontSize', fontSize*0.8, 'Color', 'k', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end

    % Drawing the force vector

    if ~isempty(forceNode)
        fx = 0;
        fy = -1;

        x0 = nodesCoordinates(forceNode, 1);
        y0 = nodesCoordinates(forceNode, 2)+1;

        quiver(x0, y0, fx, fy, 0, 'r','LineWidth', 2,'MaxHeadSize', 1);
    end

    grid on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The function below computes the static analysis (Displacement, Max Displacement and Location of Max Displacement)

function [U, U_max, Node_Max, stress, strain] = funStaticAnalysis(Load, Stiffness, ActiveDoF, N_Nodes, nodesCoordinates, connectivity, E)

    % Computing the displacement of each active DoF (others set to 0)

    U = zeros(size(Load,1),1);

    Stiffness(ActiveDoF, ActiveDoF);

    U(ActiveDoF, 1) = Stiffness(ActiveDoF,ActiveDoF)\Load(ActiveDoF,1);

    % Computing the maximal node displacement

    U_max = 0;

    for i = 1:N_Nodes

        U_Node = sqrt(U(2*i-1,1)^2 + U(2*i,1)^2);

        if U_Node > U_max

            U_max = U_Node;
            Node_Max = i;
        end
    end

    % Computing the strain and stresses of each element

    Nodes_X = nodesCoordinates(:,1);
    Nodes_Y = nodesCoordinates(:,2);

    N_Elements = size(connectivity,1);

    strain = zeros(N_Elements,1);
    stress = zeros(N_Elements,1);

    for i = 1:N_Elements

        % Extracting the i-th element parameters

        localIndices = connectivity(i,:);
    
        % Computing the initial x-length, y-length and total length

        x0 = Nodes_X(localIndices(2)) - Nodes_X(localIndices(1));
        y0 = Nodes_Y(localIndices(2)) - Nodes_Y(localIndices(1));
        L_0 = sqrt(x0*x0+y0*y0);

        % Adding the deformation

        xdef = (Nodes_X(localIndices(2)) + U(2*localIndices(2)-1)) - (Nodes_X(localIndices(1)) + U(2*localIndices(1)-1));
        ydef = (Nodes_Y(localIndices(2)) + U(2*localIndices(2))) - (Nodes_Y(localIndices(1)) + U(2*localIndices(1)));
        L_def = sqrt(xdef*xdef+ydef*ydef);

        % Computing the strain and stress

        strain(i) = (L_def-L_0)/L_0;
        stress(i) = E * strain(i);

    end
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


