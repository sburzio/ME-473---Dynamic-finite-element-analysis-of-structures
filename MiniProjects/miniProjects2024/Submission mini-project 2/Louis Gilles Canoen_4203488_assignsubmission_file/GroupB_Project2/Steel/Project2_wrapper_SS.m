%% Problem parameters

clear
close all

% Problem inputs
material = 'steel'; % options are: steel, aluminum, wood
f_base = -50000; % applied force, in N, at node 5, in y direction
scale_deformation = 2000; % multiplication factor for visualizing node deformation
nodesCoordinates_base = [0 0;-1 0;20 0; 21 0;10 0];
%nodesCoordinates_base = [4.5 0;-1 0;15.5 0; 21 0;10 0];
disp_max = 0.04; % maximum nodal displacement (m)
nf_1_range_dead = [1.5 10]; % range within which the first natural frequency must not fall (in Hz)
nf_1_goal = 50; % goal minimum first natural frequency (in Hz)


% Variables to vary/test
a = 0.1:0.1:25; % vector of heights of node 5 to test
complexity = 0:10; % complexities of truss structure to test (either look at 1:10 or 0:10)

% create empty vectors to store results
results = struct;
results_tot = zeros(length(a)*length(complexity),6);

% material properties
switch material
    case 'steel' % using 304 stainless steel  - corrosion resistent
        E = 200 * 10^9; %elastic modulus in Pa
        rho = 7850; % density in kg/m^3
        Amin = 0.0015; % minimum cross-sectional area, in m^2
        Amax = 0.0050; % maximum cross-sectional area, in m^2
        vol_lim = 1; % maximum material volume, in m^3
        ys = 215 * 10^6; % yield stress, in Pa
    case 'aluminum' % 6061-T6, commonly used
        E = 69 * 10^9; %elastic modulus in Pa
        rho = 2700; % density in kg/m^3
        Amin = 0.0020; % minimum cross-sectional area, in m^2
        Amax = 0.0070; % maximum cross-sectional area, in m^2
        vol_lim = 1.7; % maximum material volume, in m^3
        ys = 110 * 10^6; % yield stress, in MPa
    case 'wood'
        E = 12 * 10^9; %elastic modulus in Pa
        rho = 600; % density in kg/m^3
        Amin = 0.0030; % minimum cross-sectional area, in m^2
        Amax = 0.0100; % maximum cross-sectional area, in m^2
        vol_lim = 2; % maximum material volume, in m^3
        ys = 30 * 10^6; % yield stress, in MPa
end


%% Plot simple variations in truss shapes for demonstration

nodesCoordinates_base(end,2) = 5; % set height to something reasonable
buffer = 0.2; % spacing of point labels from point
bar_plot = false;

for test = 0:3 % show first four types of trusses
    [nodesCoordinates,connectivity] = Howe(test,nodesCoordinates_base);
    simpleplot(nodesCoordinates,connectivity,test,buffer,bar_plot)
end

%% Run full evaluation code - cycle through parameter combinations

A = Amax; % %x set cross-sectional area
ct = 1 ; % row number for entry into results structure

for height = a % height of node 5, at which force is applied

    nodesCoordinates_base(end,2) = height;

    for j = complexity % cycle through bridge complexity
    
    [nodesCoordinates,connectivity] = Howe(j,nodesCoordinates_base); % get out nodes coordinates and connectivity table for given parameter combination

    numberOfNodes = size(nodesCoordinates,1);

    % Test case 1 for trouble shooting
    % nodesCoordinates = [nodesCoordinates_base;5 0; 15 0; 3 a-1; 17 a-1];
    % connectivity = [1 2; 1 6; 6 7;7 3; 3 4;2 8;8 5;5 9;9 4;1 8;3 9;8 6;6 5;5 7;7 9];

    % Test case 2 for trouble shooting
    % nodesCoordinates = nodesCoordinates_base;
    % connectivity = [ 2 5; 1 5; 3 5; 4 5];
    
    % create applied loads vector
    f = zeros(length(nodesCoordinates),2);
    f(10) = f_base; 
     
    % Run calculations
    run("Project2_calc_SS.m")
    
    % save results
    results(ct).a = height;
    results(ct).complexity = j;
    results(ct).stress_MPa = maxSigma_act*10^-6;
    results(ct).displacement_mm = maxDisp_act*10^3;
    results(ct).volume_m3 = Vol_act;
    results(ct).freq_1_Hz = nat_1_act;

    row = [i, j, maxSigma_act*10^-6, maxDisp_act*10^3, Vol_act, nat_1_act];
    results_tot(ct,:) = row;

    ct = ct + 1;

    end

end

%% Plot results (displacement and stress)

a_out = [results.a];
complexity_out = [results.complexity];
f_out = [results.freq_1_Hz];
disp_out = [results.displacement_mm];
stress_out = [results.stress_MPa];
volume_out = [results.volume_m3];
[X,Y] = meshgrid(a,complexity);

f_grid = nan(length(complexity), length(a));  % rows = complexity, columns = height
stress_grid = nan(length(complexity), length(a));  % rows = complexity, columns = height
disp_grid = nan(length(complexity), length(a));  % rows = complexity, columns = height
vol_grid = nan(length(complexity), length(a));  % rows = complexity, columns = height

% create arrays of final result values for each bridge
for i = 1:length(a)
    mask1 = a_out == a(i);
    for j = 1:length(complexity)
        mask2 = complexity_out == complexity(j);
        mask_tot = mask1 & mask2;
        f_grid(j,i) = f_out(mask_tot);
        disp_grid(j,i) = disp_out(mask_tot);
        vol_grid(j,i) = volume_out(mask_tot);
        stress_grid(j,i) = stress_out(mask_tot);

    end
end

% Find where maximum frequency occurs
goal_f = f_grid>=nf_1_goal;
goal_s = stress_grid<=ys*10^-6;
best = max(max(f_grid)); 
best_idx = find(f_grid==best);
[x_star, y_star] = ind2sub(size(f_grid), best_idx);

% Find where minimum stress occurs
best_stress = min(min(stress_grid)); 
best_idx_stress = find(stress_grid==best_stress);
[x_star_s, y_star_s] = ind2sub(size(stress_grid), best_idx_stress);

if length(a)*length(complexity) > 1 % if you tested multiple cases, plot each

    % plot natural frequency contours

    figure()
    surf(X,Y,f_grid)
    hold on
    plot3(X(x_star,y_star), Y(x_star,y_star), f_grid(x_star,y_star),"pentagram", 'MarkerSize', 20, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green');
    hold off
    axis xy; colorbar;
    xlabel('Bridge height (m)','interpreter','Latex'); 
    ylabel('Truss complexity','interpreter','Latex');
    %title('First Natural frequency (Hz)','interpreter','Latex');
    shading interp
    colormap('jet')
    h = colorbar;
    h.Title.String = "Hz";
    fontsize(16,"points")
    
    figure()
    contourf(X,Y,f_grid)
    hold on
    plot(X(x_star,y_star), Y(x_star,y_star), 'pentagram', 'MarkerSize', 20, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green');
    hold off
    colormap('hot'); % or 'jet', 'hot', 'cool', etc.
    xlabel('Bridge height (m)','interpreter','Latex')
    ylabel('Truss complexity','interpreter','Latex')
    %title('First Natural Frequency','interpreter','Latex')
    h = colorbar;
    colormap('jet')
    h.Title.String = "Hz";
    fontsize(16,"points")

    % plot stress contours

    figure()
    surf(X,Y,stress_grid)
    hold on
    plot3(X(x_star_s,y_star_s), Y(x_star_s,y_star_s), stress_grid(x_star_s,y_star_s),"pentagram", 'MarkerSize', 20, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green');
    hold off
    axis xy; colorbar;
    xlabel('Bridge height (m)','interpreter','Latex');
    ylabel('Truss complexity','interpreter','Latex');
    %title('Max Stress (MPa)','interpreter','Latex');
    shading interp
    colormap('jet')
    xlim([0.1, 1])
    h = colorbar;
    h.Title.String = "MPa";
    fontsize(16,"points")
    
    figure()
    contourf(X,Y,stress_grid)
    hold on
    plot(X(x_star_s,y_star_s), Y(x_star,y_star_s), 'pentagram', 'MarkerSize', 20, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green');
    hold off
    colormap('hot'); % or 'jet', 'hot', 'cool', etc.
    xlabel('Bridge height (m)','interpreter','Latex')
    ylabel('Truss complexity','interpreter','Latex')
    %title('Max Stress (MPa)','interpreter','Latex')
    xlim([0.1, 0.3])
    h = colorbar;
    colormap('jet')
    h.Title.String = "MPa";
    fontsize(16,"points")

    figure()
    semilogy(a,stress_grid(2,:),'Linewidth',2)
    xlabel('Bridge height (m)','interpreter','Latex')
    ylabel('Maximum Stress (MPa)','interpreter','Latex')
    legend('Bridge Complexity: 1','interpreter','Latex')
    grid on
    fontsize(16,"points")

    % plot max displacement contours

    figure()
    surf(X,Y,disp_grid)
    hold on
    plot3(X(x_star,y_star), Y(x_star,y_star), stress_grid(x_star,y_star),"pentagram", 'MarkerSize', 20, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green');
    hold off
    axis xy; colorbar;
    xlabel('Bridge height (m)','interpreter','Latex');
    ylabel('Truss complexity','interpreter','Latex');
    %title('Max Displacement (mm)','interpreter','Latex');
    shading interp
    colormap('jet')
    %xlim([0.1, 1])
    h = colorbar;
    h.Title.String = "mm";
    fontsize(16,"points")

    figure()
    contourf(X,Y,disp_grid)
    hold on
    plot(X(x_star,y_star), Y(x_star,y_star), 'pentagram', 'MarkerSize', 20, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green');
    hold off
    colormap('hot'); % or 'jet', 'hot', 'cool', etc.
    xlabel('Bridge height (m)','interpreter','Latex')
    ylabel('Truss complexity','interpreter','Latex')
    %title('Max Displacement (mm)','interpreter','Latex')
    %xlim([0.1, 0.3])
    h = colorbar;
    colormap('jet')
    h.Title.String = "mm";
    fontsize(16,"points")

    % plot volume contours

    figure()
    surf(X,Y,vol_grid)
    hold on
    plot3(X(x_star,y_star), Y(x_star,y_star), stress_grid(x_star,y_star),"pentagram", 'MarkerSize', 20, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green');
    hold off
    axis xy; colorbar;
    xlabel('Bridge height (m)','interpreter','Latex');
    ylabel('Truss complexity','interpreter','Latex');
    %title('Volume (m^3)','interpreter','Latex');
    shading interp
    colormap('jet')
    %xlim([0.1, 1])
    h = colorbar;
    h.Title.String = "m^3";
    fontsize(16,"points")


    % best results height and stress
    figure;
    hold on;
    contourf(X, Y, vol_grid);
    greenstar = plot(X(x_star,y_star), Y(x_star,y_star), 'p', 'MarkerSize', 20, 'MarkerEdgeColor', 'g', ...
         'MarkerFaceColor', 'g');
    grid on;
    % Labels and Title
    xlabel('Height [m]', 'Interpreter', 'latex');
    ylabel('Complexity [m]', 'Interpreter', 'latex');
    h = colorbar;
    h.Title.String = "m^3";
    fontsize(16,"points")


end

%% Find optimal solution (based on maximum natural frequency)

[max_freq_1, idx] = max([results.freq_1_Hz]);

height = results(idx).a; %height of node 5, at which force is applied
nodesCoordinates_base(end,2) = height;
j = results(idx).complexity;
[nodesCoordinates,connectivity] = Howe(j,nodesCoordinates_base);
numberOfNodes = size(nodesCoordinates,1);
    
f = zeros(length(nodesCoordinates),2);
f(10) = f_base;
      
% Run calculations
run("Project2_calc_SS.m")
    
% Perform checks
    
if maxSigma_act > ys
    fprintf('\n Yield stress failure')
elseif maxDisp_act > disp_max
    fprintf('\n Displacement failure')
elseif Vol_act > vol_lim
    fprintf('\n Volume failure')
elseif nat_1_act >= nf_1_range_dead(1) && nat_1_act <= nf_1_range_dead(2)
    fprintf('\n First natural frequency failure')
elseif nat_1_act < nf_1_goal
    fprintf('\n Not optimal first natural frequency')
end

fprintf('\n Max stress is: %.2f MPa',maxSigma_act*10^-6)
fprintf('\n Max displacement is: %.2f mm',maxDisp_act*10^3)
fprintf('\n Volume of structure: %.2f m^3',Vol_act)
fprintf('\n First natural frequency: %.2f Hz', nat_1_act)

% Plot results
colormapplot(nodesCoordinates,connectivity,q_tot_array,newNodesCoordinates,E)
ylim([-8,12])
colormapplot(nodesCoordinates,connectivity,q_tot_array,nodesCoordinates,E)
ylim([-8,12])
simpleplot(nodesCoordinates,connectivity,j,buffer,true)
xlim([-1,21])


%% Plot Mode Shapes

bar_plot = false;
col = hsv(3); 
legendHandles = gobjects(size(col,1), 1);

% first plot basic nodes
figure()
n = get(gcf,'Number');
handle_leg = noNumPlot(nodesCoordinates,connectivity,'k',n);
legendHandles(1) = handle_leg;

% then plot mode shapes over top
for modeNumber = 1:3

us = 1:2:2*numberOfNodes-1;
vs = 2:2:2*numberOfNodes;
full_modes = zeros(GDof,length(omega));
full_modes(activeDof,:) = modal_matrix;
XX = full_modes(us,modeNumber); 
YY = full_modes(vs,modeNumber);
dispNorm = max(sqrt(XX.^2+YY.^2));
scaleFact = 1e3*dispNorm;

simpleplot(nodesCoordinates+scaleFact*[XX YY],connectivity,j,buffer,bar_plot)
grid on
grid minor

handle_leg = noNumPlot(nodesCoordinates+scaleFact*[XX YY],connectivity,col(modeNumber,:),n);
legendHandles(modeNumber+1) = handle_leg;

end

legend_str = [{'Undeformed'},{'Mode 1'},{'Mode 2'},{'Mode 3'}];
legend(legendHandles, legend_str,'interpreter','latex');


