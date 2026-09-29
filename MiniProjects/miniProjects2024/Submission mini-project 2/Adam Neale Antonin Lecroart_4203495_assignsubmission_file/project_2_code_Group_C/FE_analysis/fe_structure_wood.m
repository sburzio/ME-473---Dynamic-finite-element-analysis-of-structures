%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FE-structure-wood-bridge                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script defines the mesh of a finite element structure to be
% analyzed for the defined wooden bridge. 

% Applied load [N]:
f0 = 50e3; 

% Bridge limitations:
allowable_volume = 2;               % Maximum allowable volume [m^3]
allowable_stress = 40e6;            % Maximum allowable stress [Pa]

% Define properties for elements:
area = [88,34,88,72,78,78,72,32,30,30,32,93,93,74,74]'*1e-4;    % [m^2]
young = ones(length(area),1)*12e9;                              % [Pa]
density = ones(length(area),1)*600;                             % [kg/m^3]

% Define property matrix:
% prop(elem,:) = [Area [m^2], Young's modulus [Pa], Density [kg/m^3]] 
prop = [area,young,density];

% Node coordinates X(node,:) = [x,y] (in meters)
X = [0,0; 
     -1,0;
     20,0;
     21,0;
     10,5;
     4,0;
     1.5,2;
     16,0;
     18.5,2
     ];

% Define connectivity matrix connectivity(elem,:)=[node1, node2]
connectivity = [1,6,;
                6,8;
                8,3;
                1,7;
                7,5;
                5,9;
                9,3;
                7,6;
                6,5;
                5,8;
                8,9;
                2,7;
                4,9;
                1,5;
                5,3];

% Define load vector [N]:
r = [zeros(1,9),-f0,zeros(1,8)]'; 

% Define constrained degrees of freedom (DOF):
constrained_dof = [1,2,3,4,5,6,7,8];





