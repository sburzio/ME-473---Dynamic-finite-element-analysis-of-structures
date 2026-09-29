%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     FE-structure-aluminium-bridge                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script defines the mesh of a finite element structure to be
% analyzed for the defined aluminium bridge. 

% Applied load [N]:
f0 = 50e3; 

% Bridge specific limitations:
allowable_volume = 1.7;                 % Maximum allowable volume [m^3]
allowable_stress = 95e6;                % Maximum allowable stress [Pa]

% Define properties for elements:
area = [57,64,23,55,23,53,22,20,25,61,22,57,64,23,55,23,53,22,20,25,...
        61]'*1e-4;                                  % [m^2]
young = ones(length(area),1)*69e9;                  % [Pa]
density = ones(length(area),1)*2700;                % [kg/m^3]

% Define property matrix:
% prop(elem,:) = [Area [m^2], Young's modulus [Pa], Density [kg/m^3]] 
prop = [area,young,density];

% Node coordinates X(node,:) = [x,y] (in meters)
X = [0,0; 
     -1,0;
     20,0;
     21,0;
     10,5;
     1.9,0;
     4,0;
     2.6,3;
     6,4.4;
     10,0;
     18.1,0;
     16,0;
     17.4,3;
     14,4.4;
     ];

% Define connectivity matrix connectivity(elem,:)=[node1, node2]
connectivity = [1,6;
                1,8;
                8,6;
                6,7;
                8,7;
                8,9;
                9,7;
                9,10;
                7,10;
                9,5;
                5,10;
                3,11;
                3,13;
                13,11;
                11,12;
                13,12;
                13,14;
                14,12;
                14,10;
                12,10;
                14,5];

% Define load vector [N]:
r = [zeros(1,9),-f0,zeros(1,18)]'; 

% Define constrained degrees of freedom (DOF):
constrained_dof = [1,2,3,4,5,6,7,8];
