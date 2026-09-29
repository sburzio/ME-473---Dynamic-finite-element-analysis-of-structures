% Mini Project 2 processing

%% Material properties

% EA: axial stiffness
EA = E*A; 
% rhoA: mass per unit length 
rhoA = rho*A;

%% Geometry and Mesh Connectivity

numberOfNodes = size(nodesCoordinates,1);
GDof = 2*numberOfNodes;

%% Global Stiffness and Mass Matrices

% compute and assemble the structure stiffness matrix
stiffness = formStiffness2Dtruss(GDof, connectivity, nodesCoordinates, EA);
% compute and assemble the structure maass matrix
mass = formConsistentMass2Dtruss(GDof,connectivity, nodesCoordinates, rhoA);

%% Applying boundary conditions and solving eigenproblem

nodes_prescribed = 1:4;
nodes_free = setdiff([1:numberOfNodes],transpose(nodes_prescribed));
prescribedDof = transpose(1:2*length(nodes_prescribed)); % nodes 1 - 4 are fixed so first 8 DoF are fixed
activeDof = setdiff(transpose((1:GDof)), prescribedDof);
stiffness_freeDofs = stiffness(activeDof,activeDof);
mass_freeDofs = mass(activeDof,activeDof);
f_freeDofs = f(activeDof);
[modal_matrix,omega] = computeFrequenciesAndModes(GDof,prescribedDof,stiffness,mass,0);
frequencies = omega./(2*pi);

%% Mode Normalization and Orthogonality Checks

for modeNumber = 1:size(modal_matrix,2)
    norm = sqrt(transpose(modal_matrix(:,modeNumber))*mass_freeDofs*modal_matrix(:,modeNumber));
end

% Check orthogonality with respect to stiffness matrix (should give eigenvalues):
orthocheck_eigen = transpose(modal_matrix)*stiffness_freeDofs*modal_matrix;
%frequencies.^2;

% Check orthogonality with respect to mass matrix (should give identity matrix):
orthocheck_identity = transpose(modal_matrix)*mass_freeDofs*modal_matrix;

%% Static analysis

% Static analysis: Solve for q directly
q_tot_array = zeros(numberOfNodes,2);
q_static = stiffness_freeDofs \ f_freeDofs;

q_array = reshape(q_static,2,length(nodes_free))'; % only non-fixed nodes
q_tot_array(nodes_free,:) = q_array; % all nodes, including fixed
newNodesCoordinates = nodesCoordinates;
newNodesCoordinates(nodes_free,:) = nodesCoordinates(nodes_free,:) + q_array*scale_deformation;

%% Get out stresses in each element

numberOfElements = size(connectivity,1);
nodeCoordinateX = nodesCoordinates(:,1);
nodeCoordinateY = nodesCoordinates(:,2);
total_element_length = 0;
sigma_tot = zeros(numberOfElements,1);

for i = 1:numberOfElements
    localIndices = connectivity(i,:);   
    xa = nodeCoordinateX(localIndices(2))-nodeCoordinateX(localIndices(1));
    ya = nodeCoordinateY(localIndices(2))-nodeCoordinateY(localIndices(1));
    elementLength = sqrt(xa*xa+ya*ya);
    C = xa/elementLength;
    S = ya/elementLength;

    total_element_length = total_element_length + elementLength;
    q_loc1 = q_tot_array(localIndices(1),:);
    q_loc2 = q_tot_array(localIndices(2),:);
    q_loc_tot = [q_loc1'; q_loc2'];

    sigma_tot(i) = E/elementLength*[-C -S C S]*q_loc_tot; % stress in Pascals

end

%% Outputs

displacement = vecnorm(q_array');
maxSigma_act = max(abs(sigma_tot));
maxDisp_act = max(displacement);
Vol_act = A*total_element_length;
nat_1_act = frequencies(1);


