% Mini-Project 1 Code: Natural frequencies and mode shapes of a uniform bar

% ME-473: Dynamic finite element analysis of structures

% Group F:
% Adrien MAITROT
% Gaël GERBER
% Marko MITRIC
% Hippolyte SOULIER

% Final version done 20th March 2025

clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters

n = 1;      % [-]       -   !!!! Number of elements !!!!

L = 1;      % [m]       -   Beam length
A = 0.01;   % [m^2]     -   Beam cross-section area
E = 210e9;  % [Pa]      -   Young's modulus
rho = 7850; % [kg/m^3]  -   Density
m = 3*n+1;  % [-]       -   Number of nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defining the symbolic parameters e which is the dimensionless coordinates
% of the master element, and lambda the eigenvalue

syms e 

% Defining the shape function matrix and its gradient

H(:,1) = (1/16)*(1-9*e^2)*(e-1);
H(:,2) = (9/16)*(1-e^2)*(1-3*e);
H(:,3) = (9/16)*(1-e^2)*(3*e+1);
H(:,4) = (1/16)*(9*e^2-1)*(e+1);
B = diff(H,e);

% Computing the global stiffness and the mass matrix

[K, M] = Stiff_Mass_Mat(E, A, rho, e, L, n, H, B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the eigenvalues in diagonal matrix D and the eigenvectors V

[V, D] = eig(M\K); % M\K = inv(M)*K but avoids singularity issues

D = sqrt(D);

Om_1to4 = mink(maxk(D,1,1),4); % Finds the 4 smallest omega

mode_shape = zeros(length(V)+2,4);

for i = 1:4
    if(n<2 && i>2) 
        mode_shape(:,3:4) = []; % if n=1, we cannot have more than 2 modes
        break   
    end
    [r, c] = find(D == Om_1to4(i)); % computes index of omega in diagonal matrix D

    mode_shape(:,i) = [0;V(:,c);0]; % adding the BCs u(0)=u(L)=0
end

fprintf('For n = %i elements, we obtain: \n', n)
fprintf('\nOmega 1 numerical = %.0f rad/s', Om_1to4(1))
fprintf('\nOmega 2 numerical = %.0f rad/s', Om_1to4(2))

if(n>1)
    fprintf('\nOmega 3 numerical = %.0f rad/s', Om_1to4(3))
    fprintf('\nOmega 4 numerical = %.0f rad/s', Om_1to4(4))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing the analytical frequencies

om(1) = pi*sqrt(E/rho*L^2);
om(2) = 2*pi*sqrt(E/rho*L^2);
om(3) = 3*pi*sqrt(E/rho*L^2);
om(4) = 4*pi*sqrt(E/rho*L^2);

fprintf('\n\n------------------- \n')
fprintf('\nOmega 1 analytic = %.0f rad/s', om(1))
fprintf('\nOmega 2 analytic = %.0f rad/s', om(2))
fprintf('\nOmega 3 analytic = %.0f rad/s', om(3))
fprintf('\nOmega 4 analytic = %.0f rad/s', om(4))

% Computing the analytical mode shape

xx = linspace(0,L,200);

mx = linspace(0,L,m); % x-coordinates of the node

an_mode1 = sin(pi*xx/L);
an_mode2 = sin(2*pi*xx/L);
an_mode3 = sin(3*pi*xx/L);
an_mode4 = sin(4*pi*xx/L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting the results analytical vs numerical

elements = int2str(n);

figure()
plot(xx,an_mode1,'k--','LineWidth',2)
hold on
title('Mode shape 1')
xlabel('Position [m]')
axis([0 L -1.5 1.5])
A1 = sin(pi*mx(2)/L)./mode_shape(2,1);
plot(mx,A1*mode_shape(:,1),'o-b','LineWidth',2)
yline(0,'k--')
legend('Analytical', elements + " elements")
str1 = "Mode1_n=" + elements + ".png";
%exportgraphics(gcf,str1,'Resolution',600)
hold off

figure()
plot(xx,an_mode2,'k--','LineWidth',2)
hold on
title('Mode shape 2')
xlabel('Position [m]')
axis([0 L -1.5 1.5])
A2 = sin(2*pi*mx(2)/L)./mode_shape(2,2);
plot(mx,A2*mode_shape(:,2),'o-b','LineWidth',2)
yline(0,'k--')
legend('Analytical', elements + " elements")
str2 = "Mode2_n=" + elements + ".png";
%exportgraphics(gcf,str2,'Resolution',600)
hold off

if n>1
    figure()
    plot(xx,an_mode3,'k--','LineWidth',2)
    hold on
    title('Mode shape 3')
    xlabel('Position [m]')
    axis([0 L -1.5 1.5])
    A3 = sin(3*pi*mx(2)/L)./mode_shape(2,3);
    plot(mx,A3*mode_shape(:,3),'o-b','LineWidth',2)
    yline(0,'k--')
    legend('Analytical', elements + " elements")
    str3 = "Mode3_n=" + elements + ".png";
    %exportgraphics(gcf,str3,'Resolution',600)
    hold off

    figure()
    plot(xx,an_mode4,'k--','LineWidth',2)
    hold on
    title('Mode shape 4')
    xlabel('Position [m]')
    axis([0 L -1.5 1.5])
    A4 = sin(4*pi*mx(2)/L)./mode_shape(2,4);
    plot(mx,A4*mode_shape(:,4),'o-b','LineWidth',2)
    yline(0,'k--')
    legend('Analytical', elements + " elements")
    str4 = "Mode4_n=" + elements + ".png";
    %exportgraphics(gcf,str4,'Resolution',600)
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error computations

n_total = 16;

Error_n = zeros(4,n_total);
Error_norm_n = zeros(1,n_total);

for ii = 1:n_total
    
    % Computing the K and M matrices

    [Ki, Mi] = Stiff_Mass_Mat(E, A, rho, e, L, ii, H, B);

    % Computing the natural frequencies and sorting the 4 smallest

    [Vi, Di] = eig(Mi\Ki); 
    Di = sqrt(Di);
    Omega_i = mink(maxk(Di,1,1),4); 

    if ii<2
        % Computing the error and its norm

        Error_n(1:2, ii) = vpa(abs((Omega_i-om(1:2))./om(1:2)),4);
        Error_norm_n(ii) = norm(Error_n(1:2, ii));
    else
        % Computing the error and its norm
        
        Error_n(:,ii) = vpa(abs((Omega_i-om)./om),4);
        Error_norm_n(ii) = norm(Error_n(:,ii));
    end
end

% Plotting the errors

figure()
semilogy(1:n_total,Error_n(1,:),'-x',1:n_total,Error_n(2,:),'-x',2:n_total,Error_n(3,2:end),'-x',2:n_total,Error_n(4,2:end),'-x') 
xlabel('Number of elements') 
ylabel('Error of each mode [-]')
title('Convergence of the error for each mode')
legend({'First mode','Second mode','Third mode','Fourth mode'},'Location','best')
grid on

figure()
semilogy(1:n_total,Error_norm_n,'r-x') 
xlabel('Number of elements')
ylabel('Norm of the error [-]')
title('Convergence for the norm of the error')
grid on

fprintf('\n\n------------------- \n')
fprintf('\nFrequencies error for n = 1, 2, 4, 8 and 16 elements:\n')
fprintf('\nMode 1 error: %.3e | %.3e | %.3e | %.3e | %.3e', Error_n(1,[1,2,4,8,16]))
fprintf('\nMode 2 error: %.3e | %.3e | %.3e | %.3e | %.3e', Error_n(2,[1,2,4,8,16]))
fprintf('\nMode 3 error: N/A       | %.3e | %.3e | %.3e | %.3e', Error_n(3,[2,4,8,16]))
fprintf('\nMode 4 error: N/A       | %.3e | %.3e | %.3e | %.3e\n', Error_n(4,[2,4,8,16]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function

% Computing and assembling the stiffness and mass matrices for the
% discretized system

function [K, M] = Stiff_Mass_Mat(E, A, rho, e, L, n, H, B)

    % Computing the element length and the number of node

    eL = L/n;  
    m = 3*n+1;  

    % Computing the element stiffness and mass matrices

    eK = int(2*E*A/eL*(B'*B),e,[-1,1]);
    eM = int(eL*rho*A/2*(H'*H),e,[-1,1]);

    % Assembly of the element matrices into the global 

    K = zeros(m,m);
    M = zeros(m,m);
    for i = 1:n
        s = 3*(i-1)+1:3*i+1;
        K(s,s) = K(s,s)+eK;
        M(s,s) = M(s,s)+eM;
    end
    K = K(2:end-1,2:end-1);
    M = M(2:end-1,2:end-1);
end