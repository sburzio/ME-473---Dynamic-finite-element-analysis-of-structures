%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Script for project 2 in EPFL course ME-473 by group C           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs a static an dynamic finite element analysis for the
% three different bridge structures attached together with this script;
% i.e. a 2D steel bridge truss structure, a 2D alumiunum bridge truss
% structure and a 2D wooden bridge truss structere. In addition to the
% FE-analysis' the script also checks whether predfined limits are overheld
% for the analysed structure, and does also plot the undeformed structure
% together with modeshapes.

% The analysed bridge type is chosen in the initialization of the script.
% Function are placed in the bottom of the script.

%-------------------------  Initialization  ------------------------------%
clear all;
close all;

% Define which bridge to analyze (1 = steel, 2 = aluminum, 3 = wood): 
bridgeType = 1;

% Define general parameters for all bridges:
n_modes = 3;                    % Examined number of modes, 0 = all modes.
allowable_freq = 50;            % Minimum allowable freq. for mode 1 [Hz]
allowable_displ = 40e-3;        % Maximum allowable deformation [m]

% Load finite element mesh for chosen bridge
disp("-------------------------------------------------------------------")
fprintf("Chosen bridge: ")
switch (bridgeType)
    case 1
        fe_structure_steel              % Loading file for steel bridge
        fprintf("Construction steel.\n")
    case 2
        fe_structure_aluminium          % Loading file for aluminium bridge
        fprintf("Aluminium bridge.\n")
    case 3
        fe_structure_wood               % Loading file for wood bridge
        fprintf("Engineering wood bridge.\n")
    otherwise
        fprintf("Invalid bridge input. Script stopped.\n")
        return
end
disp("-------------------------------------------------------------------")

% Define mesh quantities
n_elem = size(connectivity,1);          % Number of elements
n_node = size(X,1);                     % Number of nodes
gdof = n_node*2;                        % Global degrees of freedom (dof)

%--------------------------- Compute matrices ----------------------------%

% Determine global stiffness and consistent mass matrices
[K,M] = form_stiffness_and_mass_2Dtruss(X,connectivity,prop,gdof,n_elem);

% Enforce boundary conditions
dof_active = setdiff(1:gdof,constrained_dof);   % Define active dof
K_reduced = K(dof_active,dof_active);           % Reduced stiffness matrix
M_reduced = M(dof_active,dof_active);           % Reduced mass matrix
r_reduced = r(dof_active);                      % Reduced load vec

%--------------------------- Static analysis -----------------------------%

% Determine displacement (displ.):
q_reduced = K_reduced \ r_reduced;              % Reduced displ. vec
q = zeros(gdof,1);                              % Preparing displ. vec      
q(dof_active) = q_reduced;                      % Full displ. vec

% Determine maximum element stresses
stress = compute_element_stress(X,connectivity,prop,n_elem,q);

%---------------------------- Modal analysis -----------------------------%

% Determine desired number of natural frequencies and modeshapes
if n_modes > 0                               % Solve eigenvalue problem
    [modeshapes_reduced,Lambda] = eigs(K_reduced,M_reduced,n_modes, ...
        'smallestabs');       
else
    [modeshapes_reduced,Lambda] = eigs(K_reduced,M_reduced);
end

omega = sqrt(diag(Lambda));                 % Determine natural frequencies
freq_hz = omega/(2*pi);                     % Frequnecies in Hz

% Define full modeshapes 
modeshapes = zeros(gdof,length(freq_hz));
modeshapes(dof_active,:) = modeshapes_reduced;

%-------------------------- Checking conditions --------------------------%

% Checking whether volume is below allowable maximum volume:
check_volume(X,connectivity,prop,n_elem,allowable_volume)

% Checking whether displacement is within limit:
check_displ(q,allowable_displ)

% Checking whether stresses are within limit
check_stress(stress,allowable_stress)

% Checking whether the lowest freqeuncy is within limit:
check_freq(freq_hz,allowable_freq)

%------------------------------ Plotting ---------------------------------%

% Plot undeformed structure
draw2DtrussModified(X,connectivity,zeros(gdof,1),'SaveFigure',true, ...
        'FigureName',"Bridge_"+num2str(bridgeType)+"_undeformed")

% Defining scaling parameters for illustrating deformation
scaling = 250;

% Plot deformed structure
draw2DtrussModified(X,connectivity,q*scaling, ...
        'stresses',stress,'LineWidth',1.5,...
        'ShowUndeformed',true,'SaveFigure',true, ...
        'FigureName',"Bridge_"+num2str(bridgeType)+"_deformed")

% Redefining scaling parameters for illustrating mode shapes
scaling = 10;

% Plot "nmodes" mode shapes in different plots
for i = 1:n_modes
    draw2DtrussModified(X,connectivity,modeshapes(:,i)*scaling, ...
        'ShowUndeformed',true,'SaveFigure',true, ...
        'FigureName',"Bridge_"+num2str(bridgeType)+"_mode_"+num2str(i))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Functions                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M] = form_stiffness_and_mass_2Dtruss(X,connectivity, ...
    prop,gdof,n_elem)
    % Function determining stiffness and consistent mass matrix for a 
    %   2D truss structure
    % Input:
        % X(node,:) = [x,y]: Node coordinate matrix
        % connectivity(elem,:) = [node 1, node 2]: Connectivity matrix
        % prop(elem,:) = [area, young's modulus, density]:
        %       element property matrix
        % gdof: global degrees of freedom
        % n_elem: number of elements

    % Prepare arrays for global matrices:
    K = zeros(gdof);
    M = zeros(gdof);
    
    % Determine element stiffness and consistent mass matrices,
    % and add to global matrix. 
    for e=1:n_elem
        % Define index for element degrees of freedom
        edof = [connectivity(e,1)*2-1,connectivity(e,1)*2, ...
            connectivity(e,2)*2-1,connectivity(e,2)*2]; 

        % Define element properties:
        A = prop(e,1);         % Cross sectional area of element.
        E = prop(e,2);         % Young's modulus of element.
        rho = prop(e,3);       % Density of element.
        
        % Define length of undeformed element:
        dx = X(connectivity(e,2),1)-X(connectivity(e,1),1);   % x_2-x_1.
        dy = X(connectivity(e,2),2)-X(connectivity(e,1),2);   % y_2-y_1.
        l0 = sqrt(dx^2+dy^2);               % Length of undeformed element.
    
        % Define properties for transformation matrix:
        c  = dx/l0;                 % cos(theta)
        s  = dy/l0;                 % sin(theta)
        T  = [c^2,s*c;s*c,s^2];     % Basis part of transformation matrix.
        
        % Define element stiffness and consistent mass matrix: 
        Ke = E*A/l0*[T,-T;-T,T];        % Element stiffness matrix.
        Me = rho*A*l0/6*[2*T,T;T,2*T];  % Element consistent mass matrix.
    
        % Add element matrices to global matrices.
        K(edof,edof) = K(edof,edof) + Ke;
        M(edof,edof) = M(edof,edof) + Me;
    end
end

function [stress] = compute_element_stress(X,connectivity,prop,n_elem,q)
    % Function determining element stresses for a 2D truss structure
    % Input:
        % X(node,:) = [x,y]: Node coordinate matrix
        % connectivity(elem,:) = [node 1, node 2]: Connectivity matrix
        % prop(elem,:) = [area, young's modulus, density]:
        %       element property matrix.
        % n_elem: number of elements
        % q: Displacement vector
    
    % Prepare array for element stresses:
    stress = zeros(n_elem,1);

    % Determine element stresses
    for e=1:n_elem
        % Define index for element degrees of freedom
        edof = [connectivity(e,1)*2-1,connectivity(e,1)*2, ...
            connectivity(e,2)*2-1,connectivity(e,2)*2]; 

        % Define relevant element properties:
        E = prop(e,2);                      % Young's modulus of element.
        
        % Define length of deformed element:
        dx = (X(connectivity(e,2),1)+q(edof(3)))...
            -(X(connectivity(e,1),1)+q(edof(1)));   % x_2e-x_1e.
        dy = (X(connectivity(e,2),2)+q(edof(4)))...
            -(X(connectivity(e,1),2)+q(edof(2)));   % y_2e-y_1e.
        le = sqrt(dx^2+dy^2);               % Length of deformed element.

        % Define transformation matrix:
        c  = dx/le;                         % cos(theta)
        s  = dy/le;                         % sin(theta)
        T  = [-c,-s,c,s];                   % Transformation matrix.

        % Determine element stress: 
        stress(e) = E/le*T*q(edof);
    end
end

function check_volume(X,connectivity,prop,n_elem,allowable_volume)
    % Function checking whether the volume of the structure is within the
    % predefined allowable limit
    % Input:
        % X(node,:) = [x,y]: Node coordinate matrix
        % connectivity(elem,:) = [node 1, node 2]: Connectivity matrix
        % prop(elem,:) = [area, young's modulus, density]:
        %       Element property matrix.
        % n_elem: number of elements
        % allowable_volume: Predefined maximum allowable volume

    % Determine total volume
    V_total = determine_total_volume(X,connectivity,prop,n_elem); 
                           
    % Check whether total volume is within limit and print result. 
    if V_total>allowable_volume
        fprintf(['Total volume is above allowable volume. Total volume:'...
            ' %.3f m^3 > %.3f m^3 \n'], V_total,allowable_volume);
    else
        fprintf(['Total volume is within limit. Total volume =' ...
            ' %.3f m^3. \n'], V_total);
    end
end

function check_displ(q,allowable_displ)
    % Function checking whether the maximum displacement of the structure 
    % is within the predefined allowable limit.
    % Input:
        % q: Displacement vector
        % allowable_displ: Predefined maximum allowable displacement

    % Prepare arrays
    displ_nodes = zeros((length(q)/2),1);       % Total displ. of nodes

    % Determine total displacement of each node:
    for n = 1:length(displ_nodes)
        displ_nodes(n) = sqrt(q(n*2-1)^2+q(n*2)^2);
    end
    
    % Define index of nodes with too large displacement
    index_too_large_displ = displ_nodes > allowable_displ;

    % Check displacement cnodition and print result
    if any(index_too_large_displ~= 0)
        fprintf(['Too large displacement at the following degrees of' ...
            ' freedom: \n']);
        disp(find(index_too_large_displ==1)')
    else
        fprintf(['Displacement is within limit for all nodes. Max dis' ...
            'placement = %.2f mm. \n'],max(displ_nodes)*1e3);
    end
end

function check_stress(stress,allowable_stress)
    % Function checking whether the maximum stress of the structure 
    % is lower than the predefined limit (the yield stress). 
    % Input:
        % stess: element stress vector
        % allowable_displ: Predefined maximum allowable stress

    index_too_large_stress = abs(stress) > allowable_stress;
    if any(index_too_large_stress~= 0)
        fprintf('Too large stress in the following elements: \n');
        disp(find(index_too_large_stress==1)')
    else
        fprintf(['All stresses are below the yield stress. Max stress' ...
            ' = %.2f MPa. \n'],max(abs(stress))*1e-6);
    end
end

function check_freq(freq,allowable_freq)
    % Function checking whether the lowest frequenecy of the structure 
    % is higher than the predefined limit. 
    % Input:
        % freq: vector containing frequencies 
        % allowable_freq: Predefined minimum frequency limit for lowest
        %                 frequency.

    if freq(1)<allowable_freq
        fprintf(['1st frequency is under allowed limit. It is  %.3f Hz' ...
            ' < %.3f Hz \n'], freq(1),allowable_freq);
    else
        fprintf(['The lowest frequency is within limit. Determined ' ...
            'frequencies are: ', repmat('%g Hz, ', 1, ...
            numel(freq)-1), 'and %g Hz \n'], ...
            freq);
    end
end

function [V_total] = determine_total_volume(X,connectivity,prop,n_elem)
    % Function determining total volume of undeformed structure
    % Input:
        % X(node,:) = [x,y]: Node coordinate matrix
        % connectivity(elem,:) = [node 1, node 2]: Connectivity matrix
        % prop(elem,:) = [area, young's modulus, density]:
        %       element property matrix.
        % n_elem: number of elements
   
    % Prepare variable:
    V_total = 0;                            % Total volume
    
    % Determine undeformed element volume and add to global volume. 
    for e=1:n_elem
        % Define relevant element properties:
        A = prop(e,1);                   % Cross sectional area of element.
        
        % Define length of undeformed element:
        dx = X(connectivity(e,2),1)-X(connectivity(e,1),1);   % x_2-x_1.
        dy = X(connectivity(e,2),2)-X(connectivity(e,1),2);   % y_2-y_1.
        l0 = sqrt(dx^2+dy^2);               % Length of undeformed element.

        % Define element volume:
        Ve = A*l0;
    
        % Add element volume to total volume
        V_total = V_total+Ve;
    end
end

function draw2DtrussModified(nodesCoordinates, connectivity, ...
    displacements, varargin)
% draw2DtrussModified draws a deformed 2D truss structure with node
% and element labels. The deformed structure is drawn using with linear
% shape-functions. If desired a colorbar can be added for stresses. 
% 
% Usage:
%   draw2Dtruss(nodesCoordinates, connectivity,displacements)
%   draw2Dtruss(nodesCoordinates, connectivity, 'LineColor', 'r', ...
%               'Title', 'My Truss')
%
% Inputs:
%   nodesCoordinates - Nx2 matrix of node coordinates [x, y]
%   connectivity     - Mx2 matrix of node indices defining truss elements
%   displacements    - Global displacement vector (length = 2N)
%   varargin         - (Optional) name-value pairs:
%                      'LineColor' (default: 'b')
%                      'LineWidth' (default: 1.2)
%                      'MarkerSize' (default: 18)
%                      'FontSize' (default: 14)
%                      'ShowNodeNumbers' (default: true)
%                      'Title' (default: '')
%                      'ShowUndeformed' (default: false)
%                      'SaveFigure' (default: false)
%                      'FigureName' (default: "2D-Truss")
%                      'WhitespaceFactor' (default: 0.05)
%                      'stresses' (default: zeros(size(connectivity,1),1))

    % Parse input options
    p = inputParser;
    addParameter(p, 'LineColor', 'b');
    addParameter(p, 'LineWidth', 1.2);
    addParameter(p, 'MarkerSize', 20);
    addParameter(p, 'FontSize', 16)
    addParameter(p, 'ShowNodeNumbers', true);
    addParameter(p, 'Title', '');
    addParameter(p, 'ShowUndeformed', false);  
    addParameter(p, 'SaveFigure', false);
    addParameter(p, 'FigureName', "2D-Truss");
    addParameter(p, 'WhitespaceFactor', 0.04);
    addParameter(p, 'stresses', zeros(size(connectivity,1),1))
    parse(p, varargin{:});
    
    % Extract parsed options
    lineColor = p.Results.LineColor;
    lineWidth = p.Results.LineWidth;
    markerSize = p.Results.MarkerSize;
    showNodeNumbers = p.Results.ShowNodeNumbers;
    fontSize = p.Results.FontSize;
    titleText = p.Results.Title;
    showUndeformed = p.Results.ShowUndeformed;
    SaveFigure = p.Results.SaveFigure;
    FigureName = p.Results.FigureName;
    WhitespaceFactor = p.Results.WhitespaceFactor;
    stresses = p.Results.stresses;

    % Compute deformed coordinates
    u = displacements(1:2:end);
    v = displacements(2:2:end);
    deformedCoordinates = nodesCoordinates + [u, v];

    % Find min- and maximum coordinates
    allCoordinates = [nodesCoordinates;deformedCoordinates];
    x_min = min(allCoordinates(:,1)); x_max = max(allCoordinates(:,1));
    y_min = min(allCoordinates(:,2)); y_max = max(allCoordinates(:,2));
    x_diff = abs(x_max-x_min); y_diff = abs(y_max-y_min);

    % Define whitespace around 2D-structure
    whitespace = max([x_diff,y_diff])*WhitespaceFactor;

    % Define size of figure:
    sizeParameter = 1000; % [pixels]
    if x_diff >= y_diff
        figureWidth = sizeParameter;
        figureHeight = sizeParameter*...
            ((y_diff+whitespace)/(x_diff+whitespace));
    else
        figureHeight = sizeParameter;
        figureWidth = sizeParameter*...
            ((x_diff+whitespace)/(y_diff+whitespace));
    end
    
    % Define color map in case of plotting element stress'
    if any(stresses~= 0)
        cmap = jet(256);                              % Defining color map
        colorParameter = (stresses-(min(stresses))) /...
                       (max(stresses)-min(stresses)); % Normalize to [0, 1]
        colors = interp1(linspace(0,1,size(cmap, 1)),cmap,colorParameter); 
                                            % Map parameters to colormap
    end    

    % Plot figure
    figure('Position',[100,100,100+figureWidth,100+figureHeight]);
    hold on; 
    xlabel('x [m]','FontSize',fontSize); 
    ylabel('y [m]','FontSize',fontSize);
    title(titleText, 'Interpreter', 'none');
    xlim([x_min-whitespace x_max+whitespace])
    ylim([y_min-whitespace y_max+whitespace])
    set(gca,'FontSize',fontSize-2)

    daspect([1 1 1]); 
    box on
    
    % Plot undeformed nodes
    if showUndeformed
        plot(nodesCoordinates(:,1), nodesCoordinates(:,2), 'o', ...
             'MarkerSize', markerSize * 0.3, ...
             'MarkerFaceColor', 'k', ...
             'MarkerEdgeColor', 'k', ...
             'LineWidth', 1.2);
    end

    % Plot truss elements
    for e = 1:size(connectivity, 1)
        n1 = connectivity(e, 1);
        n2 = connectivity(e, 2);

        x_undeformed = [nodesCoordinates(n1, 1), nodesCoordinates(n2, 1)];
        y_undeformed = [nodesCoordinates(n1, 2), nodesCoordinates(n2, 2)];
        x_deformed = [deformedCoordinates(n1, 1), ...
            deformedCoordinates(n2, 1)];
        y_deformed = [deformedCoordinates(n1, 2), ...
            deformedCoordinates(n2, 2)];
        
        % Undeformed bar
        if showUndeformed
            plot(x_undeformed, y_undeformed,'--k', 'LineWidth', ...
                lineWidth/1.5);
        end
        
        % Deformed bar
        if any(stresses~=0)
            plot(x_deformed, y_deformed,'-','Color', colors(e, :), ...
            'LineWidth', lineWidth); 
        else
            plot(x_deformed, y_deformed, '-', 'Color', lineColor, ...
            'LineWidth', lineWidth);
        end

        % Midpoint of deformed element
        midX = mean(x_deformed);
        midY = mean(y_deformed);

        % Plot white triangle with black edge at midpoint
        plot(midX, midY, 's', ...
            'MarkerSize', markerSize * 1.2, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'w');

        % Label the truss element number (larger, black font)
        text(midX, midY, sprintf('%d', e), ...
            'Color', 'b', ...nodes
            'FontSize', fontSize, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end

    % Plot nodes with white interior circles
    plot(deformedCoordinates(:,1), deformedCoordinates(:,2), 'o', ...
         'MarkerSize', markerSize, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

    % Node numbers inside the circle
    if showNodeNumbers
        for e = 1:size(deformedCoordinates, 1)
            text(deformedCoordinates(e,1), deformedCoordinates(e,2), ...
                sprintf('%d', e), ...
                'FontSize', fontSize, 'Color', 'k', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end
    
    % Show color bar in case of plotting element stress'
    if any(stresses~=0)
        colormap(cmap);
        clim([min(stresses) max(stresses)]*1e-6); 
                                    % Scale colorbar to stress values
        cb = colorbar;              % Display colorbar
        cb.Label.String = 'Stress [MPa]';
        cb.Label.FontSize = fontSize;
    end

    % Add grid
    grid on;
    grid minor;

    % Save the figure as a PNG-file
    if SaveFigure
        saveas(gcf,FigureName+".png");
    end
end