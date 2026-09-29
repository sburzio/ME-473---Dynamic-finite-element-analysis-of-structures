%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Script for mini-project 1  
%                                Group C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% Define element in meshes, and prepare loop arrays
meshes = [1, 2, 4, 8, 16];  % Elements in individual meshes 
n_freq = 5;                 % Number of saved frequencies
freq_cell = {};             % Cell array for saving frequencies 

% Define material and geometric properties
l = 1;                      % Length (meters)
E = 2.1e11;                 % Young's modulus (Pascals)
rho = 7850;                 % Density (kg/m^3)
A = 0.01;                   % Area (m^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Loop determining eigenfrequencies and -modes for different meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(meshes)

    
    % Define number of nodes and elements:
    ne = meshes(i);             % Number of elements
    npe = 3;                    % Number of nodes pr. element
    nn = ne*(npe-1)+1;          % Number of nodes

    % Node coordinates  
    X = (0:nn-1)*(l/(nn-1));

    % Topology matrix IX(node1,node3,node2)
    IX = cell(1, ne); 
    for e = 1:ne
      IX{e} = [1+(e-1)*2,3+(e-1)*2,2+(e-1)*2];
    end

    % Local shape functions
    syms xi
    ha(1) = xi*(xi-1)/2;
    ha(2) = xi*(xi+1)/2;
    ha(3) = 1-xi^2;
    
    % Compute Jacobian
    transf = cell(1, ne); % Create a cell array to store transformations
        % Create a cell array to store jacobian matrices
    jacobian_mat = cell(1, ne); 
        % Create a cell array to store inverses of jacobian matrices
    jacobian_inv = cell(1, ne); 
        % Create a cell array to store determinants of jacobian matrices
    jacobian_det = cell(1, ne);
    for e = 1:ne
        local_nodes_indices = IX{e};
        transf{e} = simplify(ha * X(local_nodes_indices)');
        jacobian_mat{e} = jacobian(transf{e}, xi);
        jacobian_det{e} = det(jacobian_mat{e});
        jacobian_inv{e} = inv(jacobian_mat{e});
    end
    
    % Elementary deformation, stiffness and mass matrices 
    B_elem = cell(1, ne);
    K_elem_exact = cell(1, ne);
    M_elem_exact = cell(1, ne);
    K_elem = cell(1, ne);
    M_elem = cell(1, ne);
    Ha = [];

    % Construct the archetypal shape functions matrix (1 x 4)
    for a=1:npe
        Ha = [Ha, ha(a)*eye(1)];
    end
    for e=1:ne
        % Construct the elementary deformation matrix (3 x 8)
        for a=1:npe
            B_elem{e} = [B_elem{e}, diff(ha(a)*jacobian_inv{e},xi)];
        end
        
        K_elem{e} = Gaussian_cubic(transpose(B_elem{e}) * A * E ...
            * B_elem{e} * jacobian_det{e});
        M_elem{e} = Gaussian_cubic(A* rho * transpose(Ha) * Ha ...
            * jacobian_det{e});
    end
    
    % Global stiffness and mass matrices 
    K_global = zeros(nn);
    M_global = zeros(nn);
    for e=1:ne
        indices = IX{e};
        K_global(indices,indices) = K_global(indices,indices) + K_elem{e};
        M_global(indices,indices) = M_global(indices,indices) + M_elem{e};
    end
    
        
    % Boundary conditions
      % nodes that have prescribed degree of freedom
    constrained_nodes = [1 nn];
      % prescribed displacements in x and y directions at constrained nodes 
    given_values_at_constrained_nodes = [0 0]; 
    free_nodes = setdiff((1:nn)',constrained_nodes);
    
    K_global_free_nodes = K_global(free_nodes, free_nodes);
    M_global_free_nodes = M_global(free_nodes, free_nodes);
    
    % Solve the generalized eigenvalue problem
    [eigenvectors, eigenvalues] = eig(double(K_global_free_nodes), ...
        double(M_global_free_nodes));
    
    % Extract natural frequencies (rad/s) and display the results 
    omega_approx = real(sqrt(diag(eigenvalues)));
    
    % Convert the natural frequencies (rad/s) to Hz
    frequencies_hz = omega_approx / (2 * pi);
    fprintf('Run with %.0f elements.\n',ne)
    fprintf('%.5f\n', frequencies_hz);
    fprintf('\n')

    % Save up to the first five frequencies
    if length(frequencies_hz) > 4
        freq_cell{i} = frequencies_hz(1:5).';
    else 
        freq_cell{i} = frequencies_hz(:).';
    end

    % Save the approximate mode shapes for the different meshes (up to 5)
    if i == 1
        mode_shapes_1 = [0;eigenvectors;0]/max(eigenvectors);
        mode_shapes_loc_1 = X;
    elseif i == 2
        mode_shapes_2(:,1) = [0;eigenvectors(:,1);0]/ ...
            max(abs(eigenvectors(:,1)));
        mode_shapes_2(:,2) = -[0;eigenvectors(:,2);0]/ ...
            max(abs(eigenvectors(:,2)));
        mode_shapes_2(:,3) = [0;eigenvectors(:,3);0]/ ...
            max(abs(eigenvectors(:,3)));
        mode_shapes_loc_2 = X;
    elseif i == 3
        mode_shapes_4(:,1) = [0;eigenvectors(:,1);0]/ ...
            max(abs(eigenvectors(:,1)));
        mode_shapes_4(:,2) = -[0;eigenvectors(:,2);0]/ ...
            max(abs(eigenvectors(:,2)));
        mode_shapes_4(:,3) = [0;eigenvectors(:,3);0]/ ...
            max(abs(eigenvectors(:,3)));
        mode_shapes_4(:,4) = [0;eigenvectors(:,4);0]/ ...
            max(abs(eigenvectors(:,4)));
        mode_shapes_4(:,5) = [0;eigenvectors(:,5);0]/ ...
            max(abs(eigenvectors(:,5)));
        mode_shapes_loc_4 = X;
    elseif i == 4
        mode_shapes_8(:,1) = [0;eigenvectors(:,1);0]/ ...
            max(abs(eigenvectors(:,1)));
        mode_shapes_8(:,2) = -[0;eigenvectors(:,2);0]/ ...
            max(abs(eigenvectors(:,2)));
        mode_shapes_8(:,3) = [0;eigenvectors(:,3);0]/ ...
            max(abs(eigenvectors(:,3)));
        mode_shapes_8(:,4) = [0;eigenvectors(:,4);0]/ ...
            max(abs(eigenvectors(:,4)));
        mode_shapes_8(:,5) = [0;eigenvectors(:,5);0]/ ...
            max(abs(eigenvectors(:,5)));
        mode_shapes_loc_8 = X;
    elseif i == 5
        mode_shapes_16(:,1) = [0;eigenvectors(:,1);0]/ ...
            max(abs(eigenvectors(:,1)));
        mode_shapes_16(:,2) = -[0;eigenvectors(:,2);0]/ ...
            max(abs(eigenvectors(:,2)));
        mode_shapes_16(:,3) = -[0;eigenvectors(:,3);0]/ ...
            max(abs(eigenvectors(:,3)));
        mode_shapes_16(:,4) = [0;eigenvectors(:,4);0]/ ...
            max(abs(eigenvectors(:,4)));
        mode_shapes_16(:,5) = -[0;eigenvectors(:,5);0]/ ...
            max(abs(eigenvectors(:,5)));
        mode_shapes_loc_16 = X;
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Plotting mode shapes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact mode shapes
t = 0:0.01:l;
mode_shape_1_exact = sin(pi*t/l);
mode_shape_2_exact = sin(2*pi*t/l);
mode_shape_3_exact = sin(3*pi*t/l);
mode_shape_4_exact = sin(4*pi*t/l);
mode_shape_5_exact = sin(5*pi*t/l);

% Plotting the approximate and exact mode shapes for the different meshes
% 1 element
figure('Position', [488, 242, 900, 350]);
plot(mode_shapes_loc_1,mode_shapes_1,t,mode_shape_1_exact, ...
    '--','LineWidth',1.5)
xlim([0,l])
ylim([-1.1,2.1])
yticks(0);
xticks([0,l]); % Set tick positions
xticklabels({'0', '$\ell$'}); % Custom labels
legend({'Approx. $1^{\mathrm{st}}$ mode shape', ...
    'Exact $1^{\mathrm{st}}$ mode shape'},'NumColumns', 2, ...
    'Location','north', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'ColorOrder', ...
    [0, 0.4470, 0.7410; 0, 0.4470, 0.7410],'FontSize', 11);
grid on;
saveas(gcf, 'ms1elem.png');

% 2 elements
figure('Position', [488, 242, 900, 350]);
plot(mode_shapes_loc_2,mode_shapes_2(:,1),t,mode_shape_1_exact,'--', ...
    mode_shapes_loc_2,mode_shapes_2(:,2),t,mode_shape_2_exact,'--', ...
    mode_shapes_loc_2,mode_shapes_2(:,3),t,mode_shape_3_exact,'--', ...
    'LineWidth',1.5);
ylim([-1.1,2.1])
yticks(0);
xticks([0,l]); % Set tick positions
xticklabels({'0', '$\ell$'}); % Custom labels
legend({'Approx. $1^{\mathrm{st}}$ mode shape', ...
    'Exact $1^{\mathrm{st}}$ mode shape', ...
    'Approx. $2^{\mathrm{nd}}$ mode shape', ...
    'Exact $2^{\mathrm{nd}}$ mode shape', ...
    'Approx. $3^{\mathrm{rd}}$ mode shape', ...
    'Exact $3^{\mathrm{rd}}$ mode shape'}, ...
    'NumColumns', 3,'Location','north', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', ...
    'ColorOrder', [0, 0.4470, 0.7410; 0, 0.4470, 0.7410; 0.8500, ...
    0.3250, 0.0980; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; ...
    0.9290, 0.6940, 0.1250],'FontSize', 11);
grid on;
saveas(gcf, 'ms2elem.png');

% 4 elements
figure('Position', [488, 242, 900, 350]);
plot(mode_shapes_loc_4,mode_shapes_4(:,1),t,mode_shape_1_exact,'--', ...
    mode_shapes_loc_4,mode_shapes_4(:,2),t,mode_shape_2_exact,'--', ...
    mode_shapes_loc_4,mode_shapes_4(:,3),t,mode_shape_3_exact,'--', ...
    mode_shapes_loc_4,mode_shapes_4(:,4),t,mode_shape_4_exact,'--', ...
    mode_shapes_loc_4,mode_shapes_4(:,5),t,mode_shape_5_exact,'--', ...
    'LineWidth',1.5);
ylim([-1.1,2.1])
yticks(0);
xticks([0,l]); % Set tick positions
xticklabels({'0', '$\ell$'}); % Custom labels
legend({'Approx. $1^{\mathrm{st}}$ mode shape', ...
    'Exact 1$^{\mathrm{st}}$ mode shape', ...
    'Approx. $2^{\mathrm{nd}}$ mode shape', ...
    'Exact $2^{\mathrm{nd}}$ mode shape', ...
    'Approx. $3^{\mathrm{rd}}$ mode shape', ...
    'Exact $3^{\mathrm{rd}}$ mode shape', ...
    'Approx. $4^{\mathrm{th}}$ mode shape', ...
    'Exact $4^{\mathrm{th}}$ mode shape', ...
    'Approx. $5^{\mathrm{th}}$ mode shape', ...
    'Exact $5^{\mathrm{th}}$ mode shape'}, ...
    'NumColumns', 3,'Location','north', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', ...
    'ColorOrder', [0, 0.4470, 0.7410; 0, 0.4470, 0.7410; 0.8500, ...
    0.3250, 0.0980; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; ...
    0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4940, 0.1840, ...
    0.5560; 0.4660, 0.6740, 0.1880; 0.4660, 0.6740, 0.1880],'FontSize',11);
grid on;
saveas(gcf, 'ms4elem.png');

% 8 elements
figure('Position', [488, 242, 900, 350]);
plot(mode_shapes_loc_8,mode_shapes_8(:,1),t,mode_shape_1_exact,'--', ...
    mode_shapes_loc_8,mode_shapes_8(:,2),t,mode_shape_2_exact,'--', ...
    mode_shapes_loc_8,mode_shapes_8(:,3),t,mode_shape_3_exact,'--', ...
    mode_shapes_loc_8,mode_shapes_8(:,4),t,mode_shape_4_exact,'--', ...
    mode_shapes_loc_8,mode_shapes_8(:,5),t,mode_shape_5_exact,'--', ...
    'LineWidth',1.5);
ylim([-1.1,2.1])
yticks(0);
xticks([0,l]); % Set tick positions
xticklabels({'0', '$\ell$'}); % Custom labels
legend({'Approx. $1^{\mathrm{st}}$ mode shape', ...
    'Exact 1$^{\mathrm{st}}$ mode shape', ...
    'Approx. $2^{\mathrm{nd}}$ mode shape', ...
    'Exact $2^{\mathrm{nd}}$ mode shape', ...
    'Approx. $3^{\mathrm{rd}}$ mode shape', ...
    'Exact $3^{\mathrm{rd}}$ mode shape', ...
    'Approx. $4^{\mathrm{th}}$ mode shape', ...
    'Exact $4^{\mathrm{th}}$ mode shape', ...
    'Approx. $5^{\mathrm{th}}$ mode shape', ...
    'Exact $5^{\mathrm{th}}$ mode shape'}, ...
    'NumColumns', 3,'Location','north', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', ...
    'ColorOrder', [0, 0.4470, 0.7410; 0, 0.4470, 0.7410; 0.8500, ...
    0.3250, 0.0980; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; ...
    0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4940, 0.1840, ...
    0.5560; 0.4660, 0.6740, 0.1880; 0.4660, 0.6740, 0.1880],'FontSize',11);
grid on;
saveas(gcf, 'ms8elem.png');

% 16 elements
figure('Position', [488, 242, 900, 350]);
plot(mode_shapes_loc_16,mode_shapes_16(:,1),t,mode_shape_1_exact,'--', ...
    mode_shapes_loc_16,mode_shapes_16(:,2),t,mode_shape_2_exact,'--', ...
    mode_shapes_loc_16,mode_shapes_16(:,3),t,mode_shape_3_exact,'--', ...
    mode_shapes_loc_16,mode_shapes_16(:,4),t,mode_shape_4_exact,'--', ...
    mode_shapes_loc_16,mode_shapes_16(:,5),t,mode_shape_5_exact,'--', ...
    'LineWidth',1.5);
ylim([-1.1,2.1])
yticks(0);
xticks([0,l]); % Set tick positions
xticklabels({'0', '$\ell$'}); % Custom labels
legend({'Approx. $1^{\mathrm{st}}$ mode shape', ...
    'Exact 1$^{\mathrm{st}}$ mode shape', ...
    'Approx. $2^{\mathrm{nd}}$ mode shape', ...
    'Exact $2^{\mathrm{nd}}$ mode shape', ...
    'Approx. $3^{\mathrm{rd}}$ mode shape', ...
    'Exact $3^{\mathrm{rd}}$ mode shape', ...
    'Approx. $4^{\mathrm{th}}$ mode shape', ...
    'Exact $4^{\mathrm{th}}$ mode shape', ...
    'Approx. $5^{\mathrm{th}}$ mode shape', ...
    'Exact $5^{\mathrm{th}}$ mode shape'}, ...
    'NumColumns', 3,'Location','north', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', ...
    'ColorOrder', [0, 0.4470, 0.7410; 0, 0.4470, 0.7410; 0.8500, ...
    0.3250, 0.0980; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; ...
    0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4940, 0.1840, ...
    0.5560; 0.4660, 0.6740, 0.1880; 0.4660, 0.6740, 0.1880],'FontSize',11);
grid on;
saveas(gcf, 'ms16elem.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Determining error and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate analythical frequnecies (in hertz)
freq_a_Hz = (1:n_freq)/2*sqrt(E/(rho*l^2));

% Calculate the relative error for FE-frequnecies
error_mat = [];

for i = 1:length(meshes)
    for j = 1:length(freq_cell{i})
       error_mat(i,j) = abs((freq_a_Hz(j)-freq_cell{i}(j))/(freq_a_Hz(j)));
    end
end

% Set non-existent values to NaN instead of 0. 
error_mat(error_mat==0) = NaN;
% Change unit to procent
error_mat = error_mat*100;      

% Define properties for plotting errors
lWidth = 1;
mSize = 14;
fSize1 = 16;
fSize2 = 14;

% Plotting errors
figure(Position = [100, 100, 900, 400]);
semilogy(meshes,error_mat(:,1),'.-',LineWidth=lWidth,MarkerSize=mSize)
hold on;

for i = 2:length(meshes)
    semilogy(meshes,error_mat(:,i),'.-',LineWidth=lWidth,MarkerSize=mSize)
end
hold off;
grid on;

xlabel('Number of elements $N_e$', Interpreter='latex', FontSize=fSize1);
ylabel('Relative error $\varepsilon_j$ [\%]', Interpreter='latex', ...
    FontSize=fSize1);

xticks(1:16)
yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100])
xlim([0,17])
ylim([1e-4,100])

% Define legend entries
legend_entries = {'1\textsuperscript{st} frequency', ...
    '2\textsuperscript{nd} frequency', ...
    '3\textsuperscript{rd} frequency', ...
    '4\textsuperscript{th} frequency', ...
    '5\textsuperscript{th} frequency'}; 
legend(legend_entries, Interpreter='latex', FontSize = fSize2, ...
    Location = 'northeast');

set(gca, TickLabelInterpreter='latex', ... % LaTeX font for tick labels
         FontSize=fSize2, ...              % Font size for tick labels
         Box = 'on');                      % Box around the plot

 exportgraphics(gcf, 'relative_error.png', 'Resolution', 600); % Save PNG


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cubic gaussian integration function
function integral_value = Gaussian_cubic(f_sym)
    syms xi
    
    % Nodes and weights for cubic (3rd order) gaussian integration
    Gauss_nodes = [-sqrt(3/5), 0, sqrt(3/5)];
    Gauss_weights = [5/9,8/9,5/9];
    
    % Integration
    integral_value = 0;
    for i = 1:3
        integral_value = integral_value + Gauss_weights(i) ...
                         * double(subs(f_sym, xi, Gauss_nodes(i)));
    end
end