function [mass] = formLumpedMass2Dtruss(GDof, connectivity, nodesCoordinates, rhoA)
% Computes the assembled lumped mass matrix

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
    
    elementaryMass = rhoA*elementLength/2*[C*C C*S 0 0; C*S S*S 0 0; 0 0 C*C C*S; 0 0 C*S S*S];
    
    mass(elementDof,elementDof) = mass(elementDof,elementDof) + elementaryMass;
end

end