function cost = costFcn (x, GDof, connectivity, nodesCoordinates_unkonwn_symmetry, nodesCoordinates_unkonwn_axis, nodesCoordinates_fixed, prescribedDof, rho, A, E, f, symmetry, non_symmetry, axis, settings)
%cost function for GA optimization of bridge design

%% define bridge geometry
nodesCoordinates_free = [nodesCoordinates_unkonwn_symmetry; nodesCoordinates_unkonwn_axis];

for jj = 1:size(nodesCoordinates_free,1)
nodesCoordinates_free(jj, 1) = x(2*jj-1)/10;
nodesCoordinates_free(jj, 2) = x(2*jj)/10;
end

counter = 1;
for jj = 1: size(symmetry, 1)
    A(symmetry(jj,1))=x(size(nodesCoordinates_free,1)*2+counter);
    A(symmetry(jj,2))=x(size(nodesCoordinates_free,1)*2+counter);
    counter = counter + 1;
end

for jj = 1:size(non_symmetry,1)
    A(non_symmetry(jj)) = x(size(nodesCoordinates_free,1)*2+counter);
    counter = counter+1;
end
% A = x(size(nodesCoordinates_free,1)*2+1:end);

nodesCoordinates_unkonwn_symmetry_2 = zeros(size(nodesCoordinates_unkonwn_symmetry,1),2);
for jj = 1: size(nodesCoordinates_unkonwn_symmetry,1)
        nodesCoordinates_unkonwn_symmetry_2(jj,1) = axis+(axis-nodesCoordinates_free(jj,1));
        nodesCoordinates_unkonwn_symmetry_2(jj,2) = nodesCoordinates_free(jj,2);
end

nodesCoordinates_free_sym = [nodesCoordinates_free; nodesCoordinates_unkonwn_symmetry_2];
nodesCoordinates = [nodesCoordinates_fixed; nodesCoordinates_free_sym];

%% computational part
%evaluate total volume of the bridge
volume = TotalVolume(connectivity, nodesCoordinates, A);


% compute and assemble the structure stiffness matrix
stiffness = formStiffness2Dtrussv2(GDof, connectivity, nodesCoordinates, E, A);

% compute and assemble the structure maass matrix
consistentMass = formConsistentMass2Dtrussv2(GDof,connectivity, nodesCoordinates, rho, A);

% evaluate eigenvalues and frequencies
[~,omega_consisMass] = computeFrequenciesAndModes(GDof,prescribedDof,stiffness,consistentMass,0);

freq_consis = omega_consisMass./(2*pi);


activeDof = setdiff(transpose((1:GDof)), prescribedDof);

% evaluate displacement
q = stiffness(activeDof, activeDof)\f(activeDof);

qq = zeros(size(prescribedDof,1)/2,1);
for jj = 1:size(activeDof,1)/2
    qq(5+jj,1) = sqrt(q(2*jj-1,1)^2+q(2*jj,1)^2);
end

maxq = max(qq*1000);


%evaluate stresses
stress = stresses(connectivity, nodesCoordinates, GDof, prescribedDof, q, E);

max_stress = max(norm(stress));

%% cost function

%evaluation of the cost function
cost = volume/settings.max_vol*settings.wV+mean(qq)/settings.lim_q*settings.wq-min(freq_consis)/settings.lim_freq*settings.wfreq-mean(stress)/settings.yield*settings.wstress;

% peanlty for set boundaries from the problem
if min(freq_consis)<= settings.lim_freq
    cost = cost + exp(settings.lim_freq/min(freq_consis))*settings.wfreq;
end

if maxq > settings.lim_q
    cost = cost + exp(maxq/settings.lim_q)*settings.wq;
end

if max_stress > settings.yield
    cost = cost + exp(max_stress/settings.yield)*settings.wstress;
end

if volume > settings.max_vol
    cost = cost + exp(volume/settings.max_vol)*settings.wV;
end

% if you fucked up, you'll find out here
if ~isreal(cost)
    fprintf("Cost is imaginary, taking only real part \n CostFcn = %.3f", cost)
    cost = real(cost);
end



if settings.print
    %plot final design

    figure
    draw2Dtruss(nodesCoordinates, connectivity, ...
    'LineColor', 'r', ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'ShowNodeNumbers', true);
    title("Final shape")

    % disp all decisional variable
    fprintf("Max displacement %d [mm] \n", maxq)
    fprintf("volume %d [m^3] \n", volume)
    fprintf("Max stress %d MPa \n", max_stress/1e6)
    fprintf("First modal frequency %d [Hz] \n", min(freq_consis))
    fprintf("Cost %d \n", cost) 

    fprintf("nodes coordinates \n")
    disp(nodesCoordinates)
    fprintf("areas \n")
    disp(A)
    fprintf("connectivity \n")
    disp(connectivity)
    
end

end