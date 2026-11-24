%% Threshold Jacobi rotations for a symmetric stiffness matrix
% K.-J. Bathe, Finite Element Procedures, Sec. 11.3.1, Example 11.9.
% Run this script to see every Jacobi rotation and the sorted eigenpairs.
K = [ 5 -4  1  0;
     -4  6 -4  1;
      1 -4  6 -4;
      0  1 -4  5];

tolerance = 1e-12;
[P, Lambda] = jacobi_method(K, tolerance, 30);
lambda = diag(Lambda);
n = size(K,1);

%% Verify the modal pairs
disp('Eigenvalues in ascending order:');
disp(lambda);
disp('Corresponding mode shapes (columns of P):');
disp(P);

relativeResidual = norm(K*P - P*Lambda,'fro') / norm(K,'fro');
orthogonalityError = norm(P.'*P - eye(n),'fro');
referenceError = norm(lambda - sort(eig(K)),inf);
fprintf('Relative modal residual: %.3e\n', relativeResidual);
fprintf('Orthogonality error:     %.3e\n', orthogonalityError);
fprintf('Maximum eigenvalue difference from eig(K): %.3e\n', referenceError);
assert(relativeResidual <= 10*tolerance && orthogonalityError <= 10*tolerance);
