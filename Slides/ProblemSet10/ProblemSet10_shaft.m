%% Problem set 10, Problem 1 -- Free-free torsional shaft
% Run the complete script. K and M below are dimensionless (Kbar and Mbar).
% Omega = omega*ell*sqrt(rho/G); Lambda = Omega^2. No toolbox is required.
clearvars
xi = (0:4)'/4;
n = numel(xi);
K = zeros(n);
M = zeros(n);
for e = 1:n-1
    nodes = [e,e+1];
    he = xi(e+1)-xi(e);
    K(nodes,nodes) = K(nodes,nodes) + [1,-1;-1,1]/he;
    M(nodes,nodes) = M(nodes,nodes) + he*[2,1;1,2]/6;
end
disp('Dimensionless stiffness K and consistent mass M')
disp(K)
disp(M)
assert(norm(K*ones(n,1)) < 1e-12)

%% Five-node reference and exact continuum modes
[Phi,lambda] = orderedModes(K,M);
k = (0:n-1)';
theta = k*pi/4;
lambdaFE = 96*(1-cos(theta))./(2+cos(theta));
omega = sqrt(lambda);
omegaContinuum = k*pi;
phiContinuumNodes = cos(xi*(k'*pi));
assert(max(abs(lambda-lambdaFE)) < 1e-10)
shapeMAC = zeros(n,1);
for j = 1:n
    shapeMAC(j) = massMAC(Phi(:,j),phiContinuumNodes(:,j),M);
end
disp(table(k,lambda,omega,omegaContinuum,shapeMAC))
disp('Full FE shapes (left end = 1)')
disp(Phi./Phi(1,:))

%% Two Rayleigh-Ritz reductions
Vgood = [ones(n,1), [16;11;0;-11;-16]];
Vbad = [2,0;-1,1;0,-2;-1,1;2,0];
bases = {Vgood,Vbad};
spaceNames = {'smooth twist and rigid rotation','symmetric, zero mean'};
ritzShapes = cell(1,2);
for s = 1:2
    V = bases{s};
    Kr = V'*K*V;
    Mr = V'*M*V;
    [Y,lambdaR] = orderedModes(Kr,Mr);
    P = V*Y;
    ritzShapes{s} = P./P(1,:);
    correlations = zeros(n,2);
    for i = 1:n
        for j = 1:2
            correlations(i,j) = massMAC(Phi(:,i),P(:,j),M);
        end
    end
    [~,matchedIndex] = max(correlations,[],1);
    targetMAC = diag(correlations(1:2,:));
    fprintf('\nTrial space: %s\n',spaceNames{s})
    disp('Reduced K, reduced M, and lifted shapes (left end = 1)')
    disp(Kr)
    disp(Mr)
    disp(ritzShapes{s})
    target = [0;1];
    matchedMode = matchedIndex(:)-1; % Mode numbers start at zero.
    omegaR = sqrt(lambdaR);
    disp(table(target,lambdaR,omegaR,targetMAC,matchedMode))
    fprintf('Absolute eigenvalue error against rigid mode: %.8g\n', ...
        abs(lambdaR(1)-lambda(1)))
    fprintf('Eigenvalue error against first elastic FE mode: %.6f %%\n', ...
        100*abs(lambdaR(2)/lambda(2)-1))
    if s == 1
        fprintf('Elastic Ritz frequency error vs full FE: %.6f %%\n', ...
            100*abs(omegaR(2)/omega(2)-1))
        fprintf('Elastic Ritz frequency error vs continuum: %.6f %%\n', ...
            100*abs(omegaR(2)/pi-1))
        fprintf('One-element frequency error vs continuum: %.6f %%\n', ...
            100*abs(sqrt(12)/pi-1))
    end
end

%% Shape comparison at a common left-end amplitude
xx = linspace(0,1,201);
figure('Name','Torsional modes on the five-node mesh');
for j = 1:2
    subplot(1,2,j)
    plot(xx,cos((j-1)*pi*xx),'k-','LineWidth',1.4)
    hold on
    plot(xi,Phi(:,j)/Phi(1,j),'b--','LineWidth',1.4)
    plot(xi,ritzShapes{1}(:,j),'ro','LineWidth',1.4)
    plot(xi,ritzShapes{2}(:,j),'s-.','LineWidth',1.1)
    grid on
    xlabel('x / ell')
    ylabel('Twist / left-end twist')
    title(sprintf('Target mode k = %d',j-1))
    legend('Exact continuum','Full FE','Smooth-space Ritz nodes', ...
        'Poor-space Ritz shape','Location','best')
end

function [P,lambda] = orderedModes(K,M)
    [P,D] = eig(K,M);
    [lambda,order] = sort(diag(D));
    P = P(:,order);
    tol = 100*eps(max(1,max(abs(lambda))));
    lambda(abs(lambda)<tol) = 0;  % Roundoff in the rigid eigenvalue.
    for j = 1:size(P,2)
        P(:,j) = P(:,j)/sqrt(P(:,j)'*M*P(:,j));
    end
end

function value = massMAC(u,v,M)
    value = abs(u'*M*v)^2/((u'*M*u)*(v'*M*v));
end
