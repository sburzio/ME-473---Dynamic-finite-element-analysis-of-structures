function [P, Lambda, rotationSummary] = jacobi_method(K, tolerance, maxSweeps, showSummary)
%JACOBI_METHOD Solve K*P = P*Lambda by threshold Jacobi sweeps (M = I).
%   [P,Lambda] = jacobi_method(K) uses tolerance 1e-12 and at most 30 sweeps.
%   Each accepted rotation prints its sweep, pivot (i,j), and the full
%   K_k, M_k, J_k, P_k. P_k contains the accumulated rotations in pivot
%   order. The returned columns of P are sorted by increasing eigenvalue.
%   Set showSummary to true to also display the rotation table after the run.
%   Request the third output to keep the table in the workspace. Each row
%   describes an accepted rotation, or has NaN rotation fields for a sweep
%   with no accepted pivots. PivotBefore is the entry before it is zeroed.

if nargin < 2, tolerance = 1e-12; end
if nargin < 3, maxSweeps = 30; end
if nargin < 4, showSummary = false; end
validateattributes(K, {'numeric'}, {'2d','square','real','finite','nonempty'});
validateattributes(tolerance, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(maxSweeps, {'numeric'}, {'scalar','integer','finite','positive'});
validateattributes(showSummary, {'logical','numeric'}, {'scalar','real'});
assert(showSummary == 0 || showSummary == 1, 'showSummary must be true, false, 1, or 0.');
assert(issymmetric(K), 'K must be a real symmetric matrix.');

n = size(K,1);
Kk = double(K);
Mk = eye(n);
Pk = eye(n);
k = 0;
recordSummary = showSummary || nargout >= 3;
summaryRows = zeros(0,9);
% For positive diagonal entries this is Bathe's relative coupling test.
% A small scale floor also permits zero diagonal entries and rigid modes.
scaleFloor = eps(max(abs(Kk(:))));

for sweep = 1:maxSweeps
    stepsBeforeSweep = k;
    oldDiagonal = diag(Kk);
    threshold = 10^(-2*sweep);
    for i = 1:n-1
        for j = i+1:n
            pivot = Kk(i,j);
            scaleI = sqrt(max(abs(Kk(i,i)), scaleFloor));
            scaleJ = sqrt(max(abs(Kk(j,j)), scaleFloor));
            if abs(pivot) <= threshold*scaleI*scaleJ
                continue
            end

            tau = (Kk(i,i) - Kk(j,j))/(2*pivot);
            if tau == 0
                t = 1;
            else
                t = sign(tau)/(abs(tau) + hypot(1,tau));
            end
            c = 1/hypot(1,t);
            s = t*c;
            rotation = [c -s; s c];
            pair = [i j];
            Jk = eye(n);
            Jk(pair,pair) = rotation;

            % Apply J_k on the right and J_k.' on the left.
            Kk(:,pair) = Kk(:,pair)*rotation;
            Kk(pair,:) = rotation.'*Kk(pair,:);
            Kk(i,j) = 0;
            Kk(j,i) = 0;
            Mk(:,pair) = Mk(:,pair)*rotation;
            Mk(pair,:) = rotation.'*Mk(pair,:);
            Pk(:,pair) = Pk(:,pair)*rotation;

            k = k + 1;
            if recordSummary
                summaryRows(end+1,:) = [sweep, threshold, k, i, j, pivot, atan2d(s,c), c, s];
            end
            fprintf('\nSweep %d, step %d, pivot (%d,%d), cos = %.7f, sin = %.7f\n', ...
                    sweep, k, i, j, c, s);
            fprintf('K_%d =\n', k); disp(Kk);
            fprintf('M_%d =\n', k); disp(Mk);
            fprintf('J_%d =\n', k); disp(Jk);
            fprintf('P_%d =\n', k); disp(Pk);
        end
    end

    diagonal = diag(Kk);
    scale = max(abs(diagonal), scaleFloor);
    diagonalChange = max(abs(diagonal - oldDiagonal)./scale);
    rootScale = sqrt(scale);
    couplings = triu(abs(Kk),1)./rootScale./rootScale.';
    maxCoupling = max(couplings(:));
    fprintf('End of sweep %d, diagonal:', sweep);
    fprintf(' %12.8f', diagonal);
    fprintf('   max coupling %.2e\n', maxCoupling);

    if recordSummary && k == stepsBeforeSweep
        summaryRows(end+1,:) = [sweep, threshold, nan(1,7)];
    end
    converged = maxCoupling <= tolerance && diagonalChange <= tolerance;
    if converged
        break
    end
end

if recordSummary
    rotationSummary = array2table(summaryRows, 'VariableNames', ...
        {'Sweep','Threshold','Step','i','j','PivotBefore','ThetaDeg','CosTheta','SinTheta'});
    if showSummary
        fprintf('\nJacobi summary: %d rotations in %d sweeps\n', k, sweep);
        disp(rotationSummary);
    end
end
if ~converged
    error('jacobi_method:NoConvergence', ...
          'Jacobi iteration did not converge within %d sweeps.', maxSweeps);
end
[lambda, order] = sort(diagonal);
P = Pk(:,order);
Lambda = diag(lambda);
end
