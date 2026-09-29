function [Volume] = TotalVolume(connectivity, nodesCoordinates, A)
% Computes the assembled consistent Volume matrix

Volume=0;

numberOfElements = size(connectivity,1);

nodeCoordinateX = nodesCoordinates(:,1);
nodeCoordinateY = nodesCoordinates(:,2);

for e = 1:numberOfElements
    % elementDof: element degrees of freedom (Dof)
    localIndices = connectivity(e,:);

    xa = nodeCoordinateX(localIndices(2))-nodeCoordinateX(localIndices(1));
    ya = nodeCoordinateY(localIndices(2))-nodeCoordinateY(localIndices(1));
    elementLength = sqrt(xa*xa+ya*ya);
    if size(A,1) == 1
            elementVolume = elementLength*A;
    else
            elementVolume = elementLength*A(e);
    end
    Volume = Volume+elementVolume;
end

end