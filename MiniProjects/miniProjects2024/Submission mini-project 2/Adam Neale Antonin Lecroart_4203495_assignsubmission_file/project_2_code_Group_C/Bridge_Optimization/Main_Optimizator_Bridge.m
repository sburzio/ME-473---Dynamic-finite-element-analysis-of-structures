clear 
close all
clc

%% Initialization with material properties

% Index to choose material (1=steel, 2=aluminum, 3= wood)
j = 1;

% E: modulus of elasticity
E = [200 69 12]*1e9;

%yield stress
yield = [300 150 50]*1e6;

% max volume (from the problem)
maxV = [1 1.5 2];

% A: area of cross section (LB & UB definition)
low_bound_A = [15 20 30]*1e-4;
upper_bound_A = [50 70 100]*1e-4;

% initial area for optimizer to start
A = 40e-4;

% rho density
rho = [7850 2700 600];

E = E(j);
rho = rho(j);



%% Problem set boundaries
% load from the problem (vertical to be applied on node 5)
f0 = 50000; %[N]

% axis of symmetry
axis = 10;

% limiting quantities set by the problem
settings.max_vol = maxV(j);% max volume
settings.yield = yield(j); % yield stress
settings.lim_q = 40; % (mm) max displacement
settings.lim_freq = 50; %(Hz) lowest permitted frequency

%% Geometry definition
% the geometry is defined in order to take into consideration the symmetry
% of the bridge
%
% nodesCoordinates_fixed represent the coordinates of the nodes set by the
% problem, they DO NOT move
%
% nodesCoordinates_unkonwn_axis represent the coordinate of the nodes on
% the axis, which are NOT reflected on the axis of symmetry
%
% nodesCoordinates_unkonwn_symmetry represent the coordinate of the movable
% nodes which are to reflected on the axis of symmetry
%
% symmetry represent the pair of symmetric elements
%
% non_symmetry represent the elements which are non simmetric, so element
% on or across the axis of symmetric
%
% connectivity is the connectivity table

% input to choose the reference design (j=1 Warren, 2 Pratt, 3 Reiforced
% Warren
j = 1;


nodesCoordinates_fixed = [0 0; -1 0; 20 0; 21 0; axis 5];

if j ==1
% % Warren

nodesCoordinates_unkonwn_axis = [];
nodesCoordinates_unkonwn_symmetry = [6 0; 3.5 4];
symmetry = [1 3; 4 7; 5 6; 8 11; 9 10]; % element that are symmetric
non_symmetry = [2];
connectivity = [1 6; 6 8; 8 3; 1 7; 7 5; 5 9; 9 3; 7 6; 6 5; 5 8; 8 9];
elseif j == 2
% % Pratt
nodesCoordinates_unkonwn_axis = [axis 0];
nodesCoordinates_unkonwn_symmetry = [3 0; 6 0; 3 5; 6 5];
symmetry = [1 12; 2 13; 3 14; 4 15; 5 16; 6 17; 7 18; 8 19; 9 20; 10 21]; % element that are symmetric
non_symmetry = [11];
connectivity = [1 6; 1 8; 8 6; 6 7; 8 7; 8 9; 9 7; 9 10; 7 10; 9 5; 5 10; 3 11; 3 13; 13 11; 11 12; 13 12; 13 14; 14 12; 14 10; 12 10; 14 5];
else
% Warren reinforced bridge
nodesCoordinates_unkonwn_axis = [];
nodesCoordinates_unkonwn_symmetry = [6 0; 3.5 4];
symmetry = [1 3; 4 7; 5 6; 8 11; 9 10; 12 13; 14 15]; % element that are symmetric
non_symmetry = [2];
connectivity = [1 6; 6 8; 8 3; 1 7; 7 5; 5 9; 9 3; 7 6; 6 5; 5 8; 8 9; 2 7; 4 9; 1 5; 5 3];
end

% presecribed degree of freedom, represent the DoF for the first 4 nodes
prescribedDof = transpose([1 2 3 4 5 6 7 8 ]);




%% Construction of the whole bridge
% construction of the whole node coordinate table for simmetric nodes
for jj = 1: size(nodesCoordinates_unkonwn_symmetry,1)
        nodesCoordinates_unkonwn_symmetry_2(jj,1) = axis+(axis-nodesCoordinates_unkonwn_symmetry(jj,1));
        nodesCoordinates_unkonwn_symmetry_2(jj,2) = nodesCoordinates_unkonwn_symmetry(jj,2);
end




% nodesCoordinates_free represent the coordinate of the free nodes (for
% half bridge) before simmetric constraint to applied
nodesCoordinates_free = [nodesCoordinates_unkonwn_symmetry; nodesCoordinates_unkonwn_axis];

%nodesCoordinates_free_sym represent the coordinate of the free nodes
nodesCoordinates_free_sym = [nodesCoordinates_free; nodesCoordinates_unkonwn_symmetry_2];

% nodesCoordinates represent the coordinate of all nodes
nodesCoordinates = [nodesCoordinates_fixed; nodesCoordinates_free_sym];

% external force vector
f = zeros(size(nodesCoordinates,1)*2,1);
f(5*2) =-f0;

% construction of the Area vector taking into account symmetry
A_vect = A*ones(size(connectivity,1)+size(non_symmetry,1),1);
A = zeros(size(connectivity,1),1);
counter = 1;
for jj = 1: size(symmetry, 1)
    A(symmetry(jj,1))=A_vect(counter);
    A(symmetry(jj,2))=A_vect(counter);
    counter = counter + 1;
end

for jj = 1:size(non_symmetry,1)
    A(non_symmetry(jj)) = A_vect(counter);
    counter = counter+1;
end

numberOfNodes = size(nodesCoordinates,1);

GDof = 2*numberOfNodes;


% plot reference geometry for verification
draw2Dtruss(nodesCoordinates, connectivity, ...
    'LineColor', 'r', ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'ShowNodeNumbers', true);



%% GA settings
% MaxStallGenerations -> max number of stall generation before exiting
% FunctionTolerance -> tollerance on the cost function
% Max generation -> max number of total generation
% PopulationSize -> size of the population at each iteration

% UseParallel -> use of the parallel toolbox in matlab for quicker
% simulation. It's generally not required for this code, unflag the
% variable if toolbox not installed

options = optimoptions('ga', 'MaxStallGenerations', 5, 'FunctionTolerance', ...
    1/1000, 'MaxGenerations', 100, 'PopulationSize', 500, 'PlotFcn', ...
    {'gaplotbestindiv', 'gaplotbestf'}, 'Display', ...
    'iter', 'UseParallel', false, 'UseVectorized', false);

% nVar -> number of free variable to be optimized
nVar = size(nodesCoordinates_free,1)*2;

%IntCon -> variable to be considered as integer (the position of the nodes
% are discretize to a mesh size of 0.1 m -> 10 cm
IntCon = 1:nVar;

nVar = nVar+size(A_vect, 1);

% range -> size of the square domain into which the nodes can move 
% (range = l/2)
settings.range = 2;


% weight in the cost function 
% wq -> displacement
% wfreq -> min(freq)
% wstress -> stress
% wV -> total volume

settings.wq = 1;
settings.wfreq = 1;
settings.wstress = 1;
settings.wV = 1;




% definition of the LB and UB vectors for each variable
for jj = 1:size(nodesCoordinates_free,1)
LB(jj*2-1) = (nodesCoordinates_free(jj,1)-settings.range)*10; 
LB(jj*2) = (nodesCoordinates_free(jj,2)-settings.range)*10; 
UB(jj*2-1) = (nodesCoordinates_free(jj,1)+settings.range)*10; 
UB(jj*2) = (nodesCoordinates_free(jj,2)+settings.range)*10; 
if nodesCoordinates_free(jj,1) == axis
    LB(jj*2-1) = axis*10;
    UB(jj*2-1) = axis*10;
end
if nodesCoordinates_free(jj,2) == 0
    LB(jj*2) = 0;
    UB(jj*2) = 0;
end
end

LB = [LB low_bound_A(j)*ones(size(A_vect, 1),1)'];
UB = [UB upper_bound_A(j)*ones(size(A_vect, 1),1)'];
%% GA 

%optimization through genetic algorithm
x = LB;
settings.print = true;
costFcn (x, GDof, connectivity, nodesCoordinates_unkonwn_symmetry, nodesCoordinates_unkonwn_axis, nodesCoordinates_fixed, prescribedDof, rho, A, E, f, symmetry, non_symmetry, axis, settings);


settings.print = false;
tic
fitnessfcn  = @(x) costFcn (x, GDof, connectivity, nodesCoordinates_unkonwn_symmetry, nodesCoordinates_unkonwn_axis, nodesCoordinates_fixed, prescribedDof, rho, A, E, f, symmetry, non_symmetry, axis, settings);
[x, fval, exitflag] = ga(fitnessfcn, nVar, [], [], [], [],...
    LB, UB, [] , IntCon, options);

computationalTime = toc;

fprintf("Required time %.2f (s)\n", computationalTime)

%% evaluate best individual and plot result

settings.print = true;
costFcn (x, GDof, connectivity, nodesCoordinates_unkonwn_symmetry, nodesCoordinates_unkonwn_axis, nodesCoordinates_fixed, prescribedDof, rho, A, E, f, symmetry, non_symmetry, axis, settings);
