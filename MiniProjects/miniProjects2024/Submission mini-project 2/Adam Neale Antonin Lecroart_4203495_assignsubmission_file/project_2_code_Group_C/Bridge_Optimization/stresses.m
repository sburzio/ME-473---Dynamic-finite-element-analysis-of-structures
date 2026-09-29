function [stresses] = stresses(connectivity, nodesCoordinates, GdoF, prescribedDof, q, E)
%function to evaluate the stresses in each elements from displacement
%results

activeDof = setdiff(transpose((1:GdoF)), prescribedDof);


counter = 1;
qq = zeros(size(activeDof,1)+size(prescribedDof,1),1);

for jj = 1:size(activeDof,1)+size(prescribedDof,1)
    if ismember(jj, prescribedDof)
        qq(jj,1) = 0;
    else
        qq(jj,1) = q(counter,1);
        counter = counter + 1;
    end
end


numberOfElements = size(connectivity,1);

stresses = zeros(numberOfElements,1);

nodeCoordinateX = nodesCoordinates(:,1);
nodeCoordinateY = nodesCoordinates(:,2);

for e = 1:numberOfElements
    localIndices = connectivity(e,:);

    xa = nodeCoordinateX(localIndices(2))+qq(localIndices(2)*2-1)-nodeCoordinateX(localIndices(1))+qq(localIndices(1)*2-1);
    ya = nodeCoordinateY(localIndices(2))+qq(localIndices(2)*2)-nodeCoordinateY(localIndices(1))-qq(localIndices(1)*2);
    elementLength = sqrt(xa*xa+ya*ya);
    C = xa/elementLength;
    S = ya/elementLength;

    stresses(e) = E/elementLength*(-C*qq(localIndices(1)*2-1)-S*qq(localIndices(1)*2)+C*qq(localIndices(2)*2-1)+S*qq(localIndices(2)*2));
end

end