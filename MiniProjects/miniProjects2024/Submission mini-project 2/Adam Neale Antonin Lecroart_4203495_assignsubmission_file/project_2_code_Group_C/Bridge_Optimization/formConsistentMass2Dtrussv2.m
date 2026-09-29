function [mass] = formConsistentMass2Dtrussv2(GDof, connectivity, nodesCoordinates, rho, A)
% Computes the assembled consistent mass matrix

mass=zeros(GDof);

numberOfElements = size(connectivity,1);

nodeCoordinateX = nodesCoordinates(:,1);
nodeCoordinateY = nodesCoordinates(:,2);

for e = 1:numberOfElements
    % elementDof: element degrees of freedom (Dof)
    localIndices = connectivity(e,:);
    elementDof = [2*localIndices(1)-1 2*localIndices(1) 2*localIndices(2)-1 2*localIndices(2)] ;

    xa = nodeCoordinateX(localIndices(2))-nodeCoordinateX(localIndices(1));
    ya = nodeCoordinateY(localIndices(2))-nodeCoordinateY(localIndices(1));
    elementLength = sqrt(xa*xa+ya*ya);
    C = xa/elementLength;
    S = ya/elementLength;
    
    elementaryMass = rho*A(e)*elementLength/6*[2*C*C 2*C*S C*C C*S; 2*C*S 2*S*S C*S S*S; C*C C*S 2*C*C 2*C*S;C*S S*S 2*C*S 2*S*S];
    
    mass(elementDof,elementDof) = mass(elementDof,elementDof) + elementaryMass;
end

end