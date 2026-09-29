%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Script for project 3 in EPFL course ME-473 by group C           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This cript contains the functions and code used to solve miniproject 3. 

% Mesh and element type is chosen in the initialization of the script.
% Additionally one here also defines: the number of modes to be determined,
% which mode to plot and whether to use exact or numerical integration. 

% As agreed with Stefano Burzio, the code is not provided for all plots 
% observed in the report. However a plotting fuction for plotting a 
% selected modeshape is included.

% Functions are placed in the bottom of the script.

%-------------------------  Initialization  ------------------------------%
clear all;
close all;

% Define element type (0 = AMC element, 1 = CR element)
elem_type = 1;

% Define number of elements pr. side:
n_elem_side = [2,4,8,16];

% Define element properties vector, prop = [length [m], thickness [m], 
% E-modulus [Pa], Poisson's ratio [-], Density [Kg/m^3]].
prop = [0.5, 5e-3, 210e9, 0.3, 7850];        

% Choose exact (true) or numerical integration (false)
int_exact = false;

% Choose number of modes to obtain for the structure
n_modes = 4;

% Select mode shape to be plotted. If plotted_modeshape = 0, no modeshapes
% are plotted. 
plotted_modeshape = 1;

% Prepare arrays for saving freq. and modeshapes for each mesh-size.
save_freq = cell(length(n_elem_side),1);
save_modeshapes = cell(length(n_elem_side),1);

% ------------------------------------------------------------------------%
%                              Run FE-Analysis                            %
% ------------------------------------------------------------------------%

% Run FE-analysis for each mesh size
for n = 1:length(n_elem_side)
    % -------------------------- Define mesh ---------------------------- %
    % Define mesh based on selected properties
    [X,connectivity,n_elem,n_node,n_dof,gdof,constrained_dof] ...
        = form_mesh(n_elem_side(n),elem_type,prop);

    % ------------------------ Define matrices -------------------------- %
    % Prepare vectors and matrices 
    K = zeros(gdof,gdof);                       % Global stiffness matrix
    M = zeros(gdof,gdof);                       % Global mass matrix

    % Define element degrees of freedom for first element
    edof = zeros(4*n_dof,1);
    for i = 1:4
        for j = 1:n_dof
            edof((i-1)*n_dof+j) = connectivity{1}(i)*n_dof-(n_dof-j);
        end
    end
    
    % Determine element matrices (the same for all elements):
    [Ke,Me] = form_element_matrices(X,edof,prop,elem_type,int_exact);
    
    % Assemble stiffness and mass matrix
    for e = 1:n_elem
        % Define vector with element degrees of freedom:
        for i = 1:4
            for j = 1:n_dof
                edof((i-1)*n_dof+j) = connectivity{e}(i)*n_dof-(n_dof-j);
            end
        end
    
        % Add stiffness and mass from element to global matrices
        K(edof,edof) = K(edof,edof)+Ke;
        M(edof,edof) = M(edof,edof)+Me;
    end

    % Define reduced matrices:
    dof_active = setdiff(1:gdof,constrained_dof);% Define active dof
    K_ff = K(dof_active,dof_active);             % Reduced stiffness matrix
    M_ff = M(dof_active,dof_active);             % Reduced mass matrix

    % ------------------------- Modal analysis -------------------------- %
    % Determine desired number of natural frequencies and modeshapes
    if n_modes > 0                          % Solve eigenvalue problem
        [modeshapes_ff,Lambda] = eigs(K_ff,M_ff,n_modes, ...
            'smallestabs');       
    else
        [modeshapes_ff,Lambda] = eigs(K_ff,M_ff);
    end

    omega = sqrt(diag(Lambda));             % Determine natural frequencies
    freq_hz = omega/(2*pi);                 % Frequnecies in Hz

    % Define full modeshapes 
    modeshapes = zeros(gdof,length(freq_hz));
    modeshapes(dof_active,:) = modeshapes_ff;
    
    % Save frequencies and mode shapes:
    save_freq{n} = freq_hz;
    save_modeshapes{n} = modeshapes;
end

% --------------------------- Print Results ----------------------------- %
% Determine analythical results for comparison
D = prop(3)*prop(2)^3/(12*(1-prop(4)^2));  % Flexural regidity of plate

% Prepare arrays:
freq_analytic_thin = zeros(n_modes,n_modes);
freq_analytic_thick = zeros(n_modes,n_modes);

% Determine analytical frequnecies using the expression for both the thin
% and thick plate. 
for n = 1:n_modes
    for m = 1:n_modes
        freq_analytic_thin(n,m) = pi^2*(n^2+m^2)/prop(1)^2* ...
            sqrt(D/(prop(2)*prop(5)))/(2*pi);
        freq_analytic_thick(n,m) = computeExactFrequency(prop(1), ...
            prop(1), prop(2), prop(4), prop(3), prop(5), 5/6, n, m)/(2*pi);
    end
end

% Sort from lowest to highest. 
freq_analytic_thin_sorted = sort(freq_analytic_thin(:));
freq_analytic_thick_sorted = sort(freq_analytic_thick(:));

% Print frequencies:
fprintf('---------------------------------------------\n')
fprintf('Analytical frequencies in Hz for thin plate:\n')
fprintf('---------------------------------------------\n')
fprintf('%6g\n',freq_analytic_thin_sorted(1:n_modes))
fprintf('\n')

fprintf('---------------------------------------------\n')
fprintf('Analytical frequencies in Hz for thick plate:\n')
fprintf('---------------------------------------------\n')
fprintf('%6g\n',freq_analytic_thick_sorted(1:n_modes))
fprintf('\n')


fprintf('---------------------------------------------\n')
fprintf('FE frequencies obtained in Hz for:\n')
fprintf('---------------------------------------------\n')
for n = 1:length(n_elem_side)
    fprintf('%ix%i mesh:\n',n_elem_side(n),n_elem_side(n))
    fprintf('%6g\n',save_freq{n})
    fprintf('\n')
end

% -------------------------- Plot modeshapes ---------------------------- %
% Plot selected modeshape for all defined meshes. 
if plotted_modeshape
    for n = 1:length(n_elem_side)
        % Define modeshape to be plotted:
        modeshape = save_modeshapes{n}(:,plotted_modeshape);
        % Define resolution (number of points pr. axis in plot). 
        resolution = 36;
        % Plot modeshape. 
        plot_modeshape(modeshape,plotted_modeshape,n_elem_side(n), ...
            elem_type,prop,resolution)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Functions                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,connectivity,n_elem,n_node,n_dof,gdof,constrained_dof] ...
    = form_mesh(n_elem_side,elem_type,prop)
    % Function defining the mesh af a square plate, which simply supported
    %   on the edges. 
    % Input:
        % n_elem_side: Number of elements per side
        % elem_type: Element type (0 = AMC element, 1 = CR element)
        % prop = [length [m], thickness [m], E-modulus [Pa], 
        %         Poisson's ratio [-], Density [Kg/m^3]]: 
    % Output: 
        % X: Nodal coordinates
        % connectivity: Element connectivity matrix
        % n_elem: Number of elements
        % n_node: Number of nodes
        % n_dof: Degrees of freedom pr. node
        % gdof: Global degrees of freedom
        % constrained_dof: Constrained global degrees of freedom.

    % Define number of elements, number of nodes and number of DOF's:
    n_elem = n_elem_side^2;                  % Total number of elements
    n_node_side = n_elem_side+1;             % Number of nodes pr. side
    n_node = n_node_side^2;                  % Total number of nodes
    n_dof = 3+elem_type;                     % Degrees of freedom pr. node
    gdof = n_dof*n_node;                     % Global degrees of freedom
    
    % Define node coordinates
    X = zeros(n_node,2);
    for i = 1:n_node_side
        for j = 1:n_node_side
            X((i-1)*n_node_side+j,:) = [i-1,j-1];
        end
    end
    X = X*prop(1)/(n_elem_side);
         
    % Define connectivity matrix
    connectivity = cell(n_elem,1);
    for i = 1:n_elem_side
        for j = 1:n_elem_side
            connectivity{(i-1)*(n_elem_side)+j} = [(i-1)*n_node_side+j,...
                i*n_node_side+j,i*n_node_side+j+1,(i-1)*n_node_side+j+1];
        end
    end

    % Define constrained degrees of freedom (simply supported on all edges)
    % Constrained nodes for left and right side:
    constrained_nodes_lr = [1:n_node_side, (n_node-n_elem_side):n_node];    
    
    % Constrained nodes for top and bottom:
    constrained_nodes_tb = zeros(n_node_side*2,0);
    for i = 1:n_node_side
        constrained_nodes_tb(i*2-1) = (i-1)*n_node_side+1;
        constrained_nodes_tb(i*2)   = (i)*n_node_side;
    end
    
    % Define constrained degrees of freedom.
    constrained_dof = [[constrained_nodes_lr,constrained_nodes_tb]*...
        n_dof-(n_dof-1), constrained_nodes_lr*n_dof-(n_dof-2), ...
        constrained_nodes_tb*n_dof-(n_dof-3)];
    
    % Reduce to unique indexes only. 
    constrained_dof = unique(constrained_dof);
end

function [Ke,Me] = form_element_matrices(X,edof,prop,elem_type,int_exact)
    % Function determining the element stiffness and mass matrix for a 
    %   AMC or CR plate element
    % Input:
        % X(node,:) = [x,y]: Node coordinate matrix
        % edof = [dof_elem_1,dof_elem_2,dof_elem_3,dof_elem_4]: 
        %   element degrees of freedom
        % prop = [length [m], thickness [m], E-modulus [Pa], 
        %         Poisson's ratio [-], Density [Kg/m^3]]: 
        %   element property vector
        % elem_type: 0 if AMC elem and 1 if CR elem.
        % int_exact: true if exact integration should be used
    % Output:
        % Ke: Element stiffness matrix
        % Me: Element mass matrix

    % Degrees of freedom per node
    n_dof = 3+elem_type;

    % Define half sidelengths of element 
    a = (X(edof(n_dof*2)/n_dof,1)-X(edof(n_dof)/n_dof,1))/2;
    b = (X(edof(n_dof*4)/n_dof,2)-X(edof(n_dof)/n_dof,2))/2;

    % Define master shape functions
    syms xi1 xi2

    switch elem_type                    % Define shape funtions
        case 0
            ha = AMC_shape_functions(a,b,xi1,xi2);
        case 1
            ha = CR_shape_functions(a,b,xi1,xi2);
        otherwise
            error('Invalid element type chosen. Program stopped.')
    end
    
    % Define shape function matrix and elementary deformation matrix
    Ha = [ha(1,:),ha(2,:),ha(3,:),ha(4,:)];
    Ba = [diff(Ha,xi1,xi1)/a^2;
          diff(Ha,xi2,xi2)/b^2;
          2*diff(Ha,xi1,xi2)/(a*b)];
    
    % Define elasticity matrix
    C = prop(3)/(1-prop(4)^2)*[1,prop(4),0;prop(4),1,0;0,0,(1-prop(4))/2];

    % Determine element matrices
    if not(int_exact)
        Ke = a*b*prop(2)^3/12*gaussian_quad_7th_order(Ba'*C*Ba);
        Me = a*b*prop(2)*prop(5)*gaussian_quad_7th_order(Ha'*Ha);
    else
        Ke = a*b*prop(2)^3/12*int(int(Ba'*C*Ba,xi1,-1,1),xi2,-1,1);
        Me = a*b*prop(2)*prop(5)*int(int(Ha'*Ha,xi1,-1,1),xi2,-1,1);
    end
end

function [ha] = AMC_shape_functions(a,b,xi1,xi2)
    % Function defines master shape functions for the AMC plate element
    xii = [-1,-1;1,-1;1,1;-1,1];        % Local coordinates of node i
    ha = sym(zeros(4,3));               % Initialize ha matrix
    
    % Define shape functions
    for i = 1:4                         
        ha(i,:) = [(1+xii(i,1).*xi1).*(1+xii(i,2).*xi2).*...
                     (2+xii(i,1).*xi1+xii(i,2).*xi2-xi1.^2-xi2.^2)/8;
                   b*(1+xii(i,1).*xi1).*(xii(i,2)+xi2).*(xi2.^2-1)/8;
                   -a*(xii(i,1)+xi1).*(xi1.^2-1).*(1+xii(i,2).*xi2)/8]';    
    end
end

function [ha] = CR_shape_functions(a,b,xi1,xi2)
    % Function defines master shape functions for the CR plate element
    xii = [-1,-1;1,-1;1,1;-1,1];        % Local coordinates of node i
    ha = sym(zeros(4,4));               % Initialize ha matrix

    % Define shape functions
    f = @(xi,xii) (-xii*xi.^3+3*xii*xi+2)/4;
    g = @(xi,xii) (xi.^3+xii*xi.^2-xi-xii)/4;
    for i = 1:4                             
        ha(i,:) = [f(xi1,xii(i,1))*f(xi2,xii(i,2))
                   b*f(xi1,xii(i,1))*g(xi2,xii(i,2));
                   -a*g(xi1,xii(i,1))*f(xi2,xii(i,2));
                   a*b*g(xi1,xii(i,1))*g(xi2,xii(i,2))];
    end
end

function integral_value = gaussian_quad_7th_order(f_sym)
    % gaussian_quad_7th_order: Computes the 2D integral of a symbolic
    % function over the square [-1, 1] x [-1, 1] using 7th-order Gaussian 
    % quadrature.
    %
    % Inputs:
    %   f_sym - symbolic function of two variables (e.g., f(xi1, xi2))
    %
    % Output:
    %   integral_value - approximate value of the integral

    syms xi1 xi2

    % 7rd-order Gauss-Legendre nodes and weights for [-1, 1]
    Gauss_nodes = [-sqrt(3/7+2/7*sqrt(6/5)), -sqrt(3/7-2/7*sqrt(6/5)), ...
                    sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7+2/7*sqrt(6/5))];
    Gauss_weights = [18-sqrt(30), 18+sqrt(30),18+sqrt(30),18-sqrt(30)]/36;

    integral_value = 0;
    for i = 1:length(Gauss_nodes)
        for j = 1:length(Gauss_nodes)
            % Compute contribution for each pair of nodes
            integral_value = integral_value + ...
                Gauss_weights(i) * Gauss_weights(j)... 
                * double(subs(f_sym, [xi1 xi2], ...
                    [Gauss_nodes(i), Gauss_nodes(j)]));
        end
    end
end

function omega_exact = computeExactFrequency(length_X, length_Y, h, ...
    nu, E, rho, k, n, m)
    % This function has been given during lectures/exercises in the course. 
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
    omega_squared = (Q * P33 + 2 * P12 * P23 * P13 - P22 ...
        * P13^2 - P11 * P23^2) / (rho * h * Q);

    % Natural frequency
    omega_exact = sqrt(omega_squared);
end

function plot_modeshape(modeshape,plotted_modeshape,n_elem_side, ...
    elem_type,prop,resolution)
    % Function that plots a single mode shape for a given mesh. 
    % Input:
        % modeshape: vector containing displacement at nodes for modeshape
        % plotted_modeshape: Modeshape to be plotted
        % n_elem_side: Number of elements per side
        % elem_type: Element type (0 = AMC element, 1 = CR element)
        % prop = [length [m], thickness [m], E-modulus [Pa], 
        %         Poisson's ratio [-], Density [Kg/m^3]]: 
        % resolution: Number of points plotted along each axis. Is
        % adjusted, if the number is not a multiple of n_elem_side.

    % Define mesh
    [X,connectivity,n_elem,~,n_dof,~,~] ...
    = form_mesh(n_elem_side,elem_type,prop);

    % Define half element side lengths a and b:
    a = prop(1)/(n_elem_side*2);
    b = a;
    
    % Define the number of plotted points in each direction pr. element
    n_points = ceil(resolution/n_elem_side);
    
    % Prepare arrays:
    u3h = cell(n_elem,1);
    x_coord = cell(n_elem,1);
    y_coord = cell(n_elem,1);
    
    % Determine values for displacement function u3h(x,y):
    for e = 1:n_elem
        % Define vector with element degrees of freedom:
        edof = zeros(4*n_dof,1);
        for i = 1:4
            for j = 1:n_dof
                edof((i-1)*n_dof+j) = connectivity{e}(i)*n_dof-(n_dof-j);
            end
        end
        
        % Define x and y coordinates for each element
        x_coord{e} = linspace(X(edof(n_dof)/n_dof,1), ...
            X(edof(n_dof*2)/n_dof,1),n_points);
        y_coord{e} = linspace(X(edof(n_dof)/n_dof,2), ...
            X(edof(n_dof*4)/n_dof,2),n_points);

        % Define xi values for finding uh3{e}
        xi = linspace(-1,1,n_points);
        
        % Prepare array:
        u3h{e} = zeros(n_points,n_points);
        
        % Run loop to determine uh3(xi1,xi2)
        for i = 1:n_points
            for j = 1:n_points
                % Define shape function at given xi-values
                switch elem_type
                    case 0
                        ha = AMC_shape_functions(a,b,xi(i),xi(j));
                    case 1
                        ha = CR_shape_functions(a,b,xi(i),xi(j));
                end

                % Determine u3h:
                for k=1:4
                    u3h{e}(i,j)=u3h{e}(i,j)+ha(k,:)*modeshape(...
                        edof(k*n_dof-(n_dof-1):k*n_dof));
                end
            end    
        end
        
    end

    % Plot mode shape:
    figure;
    hold on;
    for e = 1:n_elem
        [x,y] = meshgrid(x_coord{e},y_coord{e});
        surf(x,y,u3h{e}');
    end
    hold off;

    % Define plot properties
    colormap(jet);shading interp;c = colorbar;
    title(c, 'uh3(x,y) [-]');
    set(gca,'Fontsize',12);
    grid on;

    % Define labels and title
    xlabel('x [m]'); ylabel('y [m]');zlabel('uh3(x,y) [-]');
    title(['Modeshape for mode ', num2str(plotted_modeshape),' for ', ...
        num2str(n_elem_side),'x',num2str(n_elem_side),' mesh.']);
end