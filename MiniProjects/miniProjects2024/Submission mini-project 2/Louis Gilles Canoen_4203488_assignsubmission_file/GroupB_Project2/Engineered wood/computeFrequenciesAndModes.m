function [modes,omega] = computeFrequenciesAndModes(GDof,prescribedDof,stiffness,mass,maxEigenvalues)
% function to compute eigenvalues and eigenvectors
% GDof: number of degree of freedom
% prescribedDof: bounded boundary dofs
% stiffness: stiffness matrix
% mass: mass matrix
% maxEigenvalues: maximum eigenvalues to be computed. If 0 all the
% eigenvalues are requested (suggested for structures)

activeDof = setdiff(transpose((1:GDof)), prescribedDof);

if maxEigenvalues == 0
    [modes,D] = eig(stiffness(activeDof,activeDof),mass(activeDof,activeDof));
else
    [modes,D] = eigs(stiffness(activeDof,activeDof),mass(activeDof,activeDof),maxEigenvalues,'smallestabs');
end

eigenvalues = diag(D);
omega = sqrt(eigenvalues);

end