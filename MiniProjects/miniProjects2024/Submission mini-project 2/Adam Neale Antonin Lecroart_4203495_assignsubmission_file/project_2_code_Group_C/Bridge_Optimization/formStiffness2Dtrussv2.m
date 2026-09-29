function [stiffness] = formStiffness2Dtrussv2(GDof, connectivity, nodesCoordinates, E, A)
% Computes the assembled stiffness matrix

stiffness=zeros(GDof);

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
    
    elementaryStiffness = E*A(e)./elementLength*[C*C C*S -C*C -C*S; C*S S*S -C*S -S*S;-C*C -C*S C*C C*S;-C*S -S*S C*S S*S];
    
    stiffness(elementDof,elementDof) = stiffness(elementDof,elementDof) + elementaryStiffness;
end

end