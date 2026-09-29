% MATLAB code for natural frequencies of a 1D bar using quadratic elements
% with direct integration using MATLAB's integral function.
% Clamped at x = 0, free at x = L (cantilever condition).

close all
clear
clc

nElemMax = 16;
nPlot = [1, 2, 4, 8, 16]; % Elements chosen for convergence plots
modeShapes1 = cell(length(nPlot), 1);
modeShapes2 = cell(length(nPlot), 1);
modeShapes3 = cell(length(nPlot), 1);
modeShapes4 = cell(length(nPlot), 1);
modeShapes5 = cell(length(nPlot), 1);
x_values = cell(length(nPlot), 1);

% Parameters
L = 1;        % Length of the bar
E = 210e9;    % Young's modulus (GPa)
rho = 7850;   % Density (kg/m^3)
A = 0.01;     % Cross-sectional area (m^2)

for nElem = 1:nElemMax
    nElem
    nNodes = 2 * nElem + 1; % Number of nodes (quadratic elements)
    x = linspace(0, L, nNodes); % Node coordinates

    % Connectivity matrix (Each element has 3 nodes: left, right, middle)
    elements = zeros(nElem, 3);
    for e = 1:nElem
        elements(e, :) = [2*e - 1, 2*e + 1, 2*e];
    end

    % Initialize global stiffness and mass matrices
    K = zeros(nNodes, nNodes);
    M = zeros(nNodes, nNodes);

    % Element stiffness and mass matrices computation using direct integration
    syms xi;
    N = [(xi * (xi - 1)) / 2, (xi * (xi + 1)) / 2, (1 - xi^2)];
    dN_dxi = [xi - 0.5, xi + 0.5, -2 * xi];

    for e = 1:nElem
        nodes = elements(e, :);
        x1 = x(nodes(1));
        x2 = x(nodes(2)); % Rightmost node
        le = x2 - x1; % Element length
        dxdxi = le / 2;
        dN_dx = dN_dxi / dxdxi;

        % Element stiffness matrix integral
        Ke = double(int(E * A * (dN_dx' * dN_dx) * dxdxi, xi, -1, 1));

        % Element mass matrix integral
        Me = double(int(rho * A * (N' * N) * dxdxi, xi, -1, 1));

        % Assembly into global matrices
        K(nodes, nodes) = K(nodes, nodes) + Ke;
        M(nodes, nodes) = M(nodes, nodes) + Me;
    end

    % Apply boundary conditions (Clamped at x = 0, Free at x = L)
    fixedDOF = 1; % Only the first node is clamped (x = 0)
    freeDOF = setdiff(1:nNodes, fixedDOF); % All other DOFs are free

    % Extract reduced matrices
    K_reduced = K(freeDOF, freeDOF);
    M_reduced = M(freeDOF, freeDOF);

    % Solve the generalized eigenvalue problem
    [eigVecs, eigVals] = eig(double(K_reduced), double(M_reduced));
    omega = sqrt(diag(eigVals)); % Natural frequencies (rad/s)

    % Sort frequencies and mode shapes
    [omega, idx] = sort(omega);
    eigVecs = eigVecs(:, idx);

    % Store first five frequencies
    freq1(nElem) = omega(1) / (2 * pi);
    freq2(nElem) = omega(2) / (2 * pi);
    if nElem > 1
        freq3(nElem-1) = omega(3) / (2 * pi);
        freq4(nElem-1) = omega(4) / (2 * pi);
    end
    if nElem > 2, freq5(nElem-2) = omega(5) / (2 * pi); end

    % Store mode shapes for convergence plot
    if ismember(nElem, nPlot)
        i = find(nPlot == nElem);
        x_values{i} = x;  % Store the node positions

        % First mode shape (always exists)
        modeShapes1{i} = [0; eigVecs(:,1)];

        % Second mode shape (exists for n >= 1)
        if size(eigVecs, 2) > 1
            modeShapes2{i} = [0; eigVecs(:,2)];
        end

        % Third mode shape (exists for n >= 2)
        if size(eigVecs, 2) > 2
            modeShapes3{i} = [0; eigVecs(:,3)];
        end

        % Fourth mode shape (exists for n >= 2)
        if size(eigVecs, 2) > 3
            modeShapes4{i} = [0; eigVecs(:,4)];
        end

        % Fifth mode shape (exists for n >= 4)
        if size(eigVecs, 2) > 4
            modeShapes5{i} = [0; eigVecs(:,5)];
        end
    end
end


% Exact analytical frequencies for a cantilever beam (simplified approximation)
freq_real = (2 * (1:5) - 1)/(4*L) * sqrt(E/rho);

%% Convergence Plot for First Mode Shape (Normalized)
figure;
hold on;
colors = lines(length(nPlot));
for i = 1:length(nPlot)
    modeShapes1{i} = modeShapes1{i} / max(abs(modeShapes1{i})); % Normalize
    if modeShapes1{i}(2) < 0
        modeShapes1{i} = -modeShapes1{i};
    end
    plot(x_values{i}, modeShapes1{i}, 'Color', colors(i, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('$n = %d$', nPlot(i)));
end
xlabel('Position along bar $x$ [m]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Normalized displacement [-]', 'Interpreter', 'Latex', 'FontSize', 15);
legend('show', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'southeast');
grid on;

%% Convergence Plot for Second Mode Shape (Normalized)
figure;
hold on;
for i = 1:length(nPlot)
    if ~isempty(modeShapes2{i})
        modeShapes2{i} = modeShapes2{i} / max(abs(modeShapes2{i})); % Normalize
        if modeShapes2{i}(2) < 0
            modeShapes2{i} = -modeShapes2{i};
        end
        plot(x_values{i}, modeShapes2{i}, 'Color', colors(i, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('$n = %d$', nPlot(i)));
    end
end
xlabel('Position along bar $x$ [m]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Normalized displacement [-]', 'Interpreter', 'Latex', 'FontSize', 15);
legend('show', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'Southwest');
grid on;

%% Convergence Plot for Third Mode Shape (Normalized)
figure;
hold on;
for i = 1:length(nPlot)
    if ~isempty(modeShapes3{i})
        modeShapes3{i} = modeShapes3{i} / max(abs(modeShapes3{i})); % Normalize
        if modeShapes3{i}(2) < 0
            modeShapes3{i} = -modeShapes3{i};
        end
        plot(x_values{i}, modeShapes3{i}, 'Color', colors(i, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('$n = %d$', nPlot(i)));
    end
end
xlabel('Position along bar $x$ [m]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Normalized displacement [-]', 'Interpreter', 'Latex', 'FontSize', 15);
legend('show', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'Southwest');
grid on;

%% Convergence Plot for Fourth Mode Shape (Normalized)
figure;
hold on;
for i = 1:length(nPlot)
    if ~isempty(modeShapes4{i})
        modeShapes4{i} = modeShapes4{i} / max(abs(modeShapes4{i})); % Normalize
        if modeShapes4{i}(2) < 0
            modeShapes4{i} = -modeShapes4{i};
        end
        plot(x_values{i}, modeShapes4{i}, 'Color', colors(i, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('$n = %d$', nPlot(i)));
    end
end
xlabel('Position along bar $x$ [m]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Normalized displacement [-]', 'Interpreter', 'Latex', 'FontSize', 15);
legend('show', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'Southwest');
grid on;

%% Convergence Plot for Fifth Mode Shape (Normalized)
figure;
hold on;
for i = 1:length(nPlot)
    if ~isempty(modeShapes5{i})
        modeShapes5{i} = modeShapes5{i} / max(abs(modeShapes5{i})); % Normalize
        if modeShapes5{i}(2) < 0
            modeShapes5{i} = -modeShapes5{i};
        end
        plot(x_values{i}, modeShapes5{i}, 'Color', colors(i, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('$n = %d$', nPlot(i)));
    end
end
xlabel('Position along bar $x$ [m]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Normalized displacement [-]', 'Interpreter', 'Latex', 'FontSize', 15);
legend('show', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'Southwest');
grid on;

%% Error Analysis for First Natural Frequency
err1 = (freq1 - freq_real(1)) / freq_real(1);

figure
hold on
pfreq1 = plot(linspace(1,nElemMax,nElemMax), freq1, 'Linewidth', 1.2);
pfreq_real1 = yline(freq_real(1), 'k--', 'LineWidth', 1.2);
xticks([1, 2, 4, 8, 16])
axis([0, 16, 1290, 1300])
set(gca, 'XScale', 'log');
text(1.05, freq_real(1) + 0.4, sprintf('%.2f Hz', freq_real(1)), 'FontSize', 12, 'Color', 'k')
xlabel('Number of elements $n$ [-]', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('First natural frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 15)
legend('$f^{h}_1$ Numerical first natural frequency', '$f_1$ Analytical first natural frequency', 'Interpreter', 'latex', 'FontSize', 15, 'Location', 'southeast')
grid on

figure;
hold on;
plot(linspace(1, nElemMax, nElemMax), err1, 'LineWidth', 1.2);
plot(linspace(1, nElemMax, 100), 1./(linspace(1, nElemMax, 100).^4), '--', 'LineWidth', 1.2);
set(gca, 'YScale', 'log', 'XScale', 'log');
xticks([1, 2, 4, 8, 16])
xlabel('Number of elements $n$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Relative error [-]', 'Interpreter', 'Latex', 'FontSize', 15);
legend('$\epsilon_r$', '$1/n^4$', 'Interpreter', 'Latex', 'FontSize', 18);
grid on;

%% Error Analysis for Second Natural Frequency
err2 = (freq2 - freq_real(2)) / freq_real(2);

figure
hold on
pfreq2 = plot(linspace(1,nElemMax,nElemMax), freq2, 'Linewidth', 1.2);
pfreq_real2 = yline(freq_real(2), 'k--', 'LineWidth', 1.2);
xticks([1, 2, 4, 8, 16])
axis([0, 16, 3600, 4700])
set(gca, 'XScale', 'log');
text(1.05, freq_real(2) + 40, sprintf('%.2f Hz', freq_real(2)), 'FontSize', 12, 'Color', 'k')
xlabel('Number of elements $n$ [-]', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('Second natural frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 15)
legend('$f^{h}_2$ Numerical second natural frequency', '$f_2$ Analytical second natural frequency', 'Interpreter', 'latex', 'FontSize', 15, 'Location', 'southeast')
grid on

figure;
hold on;
plot(linspace(1, nElemMax, nElemMax), err1, 'LineWidth', 1.2);
plot(linspace(1, nElemMax, 100), 1./(linspace(1, nElemMax, 100).^4), '--', 'LineWidth', 1.2);
set(gca, 'YScale', 'log', 'XScale', 'log');
xticks([1, 2, 4, 8, 16])
xlabel('Number of elements $n$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Relative error [-]', 'Interpreter', 'Latex', 'FontSize', 15);
legend('$\epsilon_r$', '$1/n^4$', 'Interpreter', 'Latex', 'FontSize', 18);
grid on;

% %% Polynomial Fit for Convergence Rate
% p = polyfit(log(linspace(2, nElemMax, nElemMax-1)), log(err1(2:end)), 1);

%% Error Analysis for Third Natural Frequency
err3 = (freq3 - freq_real(3)) / freq_real(3);

figure
hold on
pfreq3 = plot(linspace(2,nElemMax,nElemMax-1), freq3, 'Linewidth', 1.2);
pfreq_real3 = yline(freq_real(3), 'k--', 'LineWidth', 1.2);
xticks([1, 2, 4, 8, 16])
axis([0, 16, 6000, 7400])
set(gca, 'XScale', 'log');
text(2.05, freq_real(3) + 60, sprintf('%.2f Hz', freq_real(3)), 'FontSize', 12, 'Color', 'k')
xlabel('Number of elements $n$ [-]', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('Third natural frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 15)
legend('$f^{h}_3$ Numerical third natural frequency', '$f_3$ Analytical third natural frequency', 'Interpreter', 'latex', 'FontSize', 15, 'Location', 'southeast')
grid on

%% Error Analysis for Fourth Natural Frequency
err4 = (freq4 - freq_real(4)) / freq_real(4);

figure
hold on
pfreq4 = plot(linspace(2,nElemMax,nElemMax-1), freq4, 'Linewidth', 1.2);
pfreq_real4 = yline(freq_real(4), 'k--', 'LineWidth', 1.2);
xticks([1, 2, 4, 8, 16])
axis([0, 16, 8000, 12000])
set(gca, 'XScale', 'log');
text(2.05, freq_real(4) + 150, sprintf('%.2f Hz', freq_real(4)), 'FontSize', 12, 'Color', 'k')
xlabel('Number of elements $n$ [-]', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('Fourth natural frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 15)
legend('$f^{h}_4$ Numerical fourth natural frequency', '$f_4$ Analytical fourth natural frequency', 'Interpreter', 'latex', 'FontSize', 15, 'Location', 'southeast')
grid on

%% Error Analysis for Fifth Natural Frequency
err5 = (freq5 - freq_real(5)) / freq_real(5);

figure
hold on
pfreq5 = plot(linspace(3,nElemMax,nElemMax-2), freq5, 'Linewidth', 1.2);
pfreq_real5 = yline(freq_real(5), 'k--', 'LineWidth', 1.2);
xticks([1, 2, 4, 8, 16])
axis([0, 16, 10000, 15000])
set(gca, 'XScale', 'log');
text(3.05, freq_real(5) + 200, sprintf('%.2f Hz', freq_real(5)), 'FontSize', 12, 'Color', 'k')
xlabel('Number of elements $n$ [-]', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('Fifth natural frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 15)
legend('$f^{h}_5$ Numerical fifth natural frequency', '$f_5$ Analytical fifth natural frequency', 'Interpreter', 'latex', 'FontSize', 15, 'Location', 'southeast')
grid on

%% All Together log_plot

figure
hold on
plot(linspace(1, nElemMax, nElemMax), err1, 'LineWidth', 1.2);
plot(linspace(1, nElemMax, nElemMax), err2, 'LineWidth', 1.2);
plot(linspace(2, nElemMax, nElemMax-1), err3, 'LineWidth', 1.2);
plot(linspace(2, nElemMax, nElemMax-1), err4, 'LineWidth', 1.2);
plot(linspace(3, nElemMax, nElemMax-2), err5, 'LineWidth', 1.2);
plot(linspace(1, nElemMax, 100), 1./(linspace(1, nElemMax, 100).^4), 'r--', 'LineWidth', 1.2);
xticks([1, 2, 4, 8, 16])
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Number of elements $n$ [-]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('Relative error $\epsilon_r$ [-]', 'Interpreter', 'Latex', 'FontSize', 15);
legend('1st natural frequency', '2nd natural frequency', '3rd natural frequency', '4th natural frequency', '5th natural frequency', '$1/n^4$', 'Interpreter', 'latex', 'Location', 'southwest', 'Fontsize', 11)
grid on
