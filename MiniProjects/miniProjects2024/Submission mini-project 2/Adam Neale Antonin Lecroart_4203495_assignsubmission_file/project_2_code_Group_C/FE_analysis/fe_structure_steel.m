%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       FE-structure-steel-bridge                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script defines the mesh of a finite element structure to be
% analyzed for the defined steel bridge. 

% Applied load [N]:
f0 = 50e3; 

% Bridge specific limitations:
allowable_volume = 1;               % Maximum allowable volume [m^3]
allowable_stress = 200e6;           % Maximum allowable stress [Pa]

% Define properties for elements:
area = [46,20,46,47,44,44,47,16,15,15,16]'*1e-4; % [m^2]
young = ones(length(area),1)*200e9;              % [Pa]
density = ones(length(area),1)*7850;             % [kg/m^3]

% Define property matrix:
% prop(elem,:) = [Area [m^2], Young's modulus [Pa], Density [kg/m^3]] 
prop = [area,young,density];

% Node coordinates X(node,:) = [x,y] (in meters)
X = [0,0; 
     -1,0;
     20,0;
     21,0;
     10,5;
     4.2,0;
     2.3,2.1;
     15.8,0;
     17.7,2.1
     ];

% Define connectivity matrix connectivity(elem,:)=[node1, node2]
connectivity = [1,6;
                6,8;
                8,3;
                1,7;
                7,5;
                5,9;
                9,3;
                7,6;
                6,5;
                5,8;
                8,9];

% Define load vector [N]:
r = [zeros(1,9),-f0,zeros(1,8)]'; 

% Define constrained degrees of freedom (DOF):
constrained_dof = [1,2,3,4,5,6,7,8];






