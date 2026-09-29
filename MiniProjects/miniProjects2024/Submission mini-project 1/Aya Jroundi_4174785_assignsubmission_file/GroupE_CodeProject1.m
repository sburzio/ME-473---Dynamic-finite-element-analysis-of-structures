format long  % Just because we want more digits
% Initialize an array to store the relative error of the first natural frequency
errfreq=zeros(5,1);
powers=0:4;


%for i = powers+1
%    errfreq(i)=fem(2^i);
%end
% Plot the relative error against the number of elements using a log-log scale
figure;
loglog(2.^powers,errfreq)
xlabel("Number of elements")
ylabel("Relative error in first natural frequency")
grid on;

% To make graphs appear for higher frequencies, comment the above code and
% just call the function for a single number of elements

fem(16)

function errfreq1 = fem(n)   % Finite Element Method (FEM) function
    % Material and geometric properties
    E = 210*10^9;  % Young's modulus (Pa)
    A = 0.01;      % Cross-sectional area (m^2)
    rho = 7850;    % Density (kg/m^3)
    l = 1;         % Length of the beam (m)
    ep = 4;        % Number of nodes per element. This can be changed to make the code run for any order of shape functions!
    
    p=(ep-1)*n+1 % Total number of nodes
    
    % Define symbolic variable for shape function computation
    syms x
    
    % Local Cubic Nodes : Define node positions in the reference element [-1,1]
    x_loc=linspace(-1,1,ep)
    
    % Local Cubic Shape Functions : (Lagrange interpolation polynomials)
    syms H B  % H for shape functions, B for their derivatives
        
    for i=1:ep
        H(:,i)=1;
        for j=1:ep
            if i~=j
                H(:,i) = H(:,i) * (x-x_loc(j))/(x_loc(i)-x_loc(j));
            end
        end
        B(:,i)= gradient(H(:,i),x); % Compute derivatives of shape functions
    end
    
    % Plot shape functions and their derivatives
    fplot(H,[-1 1])
    fplot(B,[-1 1])
    
    % Compute local mass and stiffness matrix (element-wise mass contribution)
    e_M=int(rho*A*H'*H,x,[-1 1])
    e_K=int(E*A*B'*B,x,[-1 1])
    
    % Compute the Jacobian determinant for transformation from local to global coordinates
    e_j=l/n/2
    
    % Initialize global mass and stiffness matrices
    M=zeros(p,p); 
    K=zeros(p,p);
    % Global matrix assembly: loop over elements
    for i = 1:n
        for j = 1:ep
            for k = 1:ep
                M((i-1)*(ep-1)+j,(i-1)*(ep-1)+k)=M((i-1)*(ep-1)+j,(i-1)*(ep-1)+k)+e_M(j,k)*e_j; 
                K((i-1)*(ep-1)+j,(i-1)*(ep-1)+k)=K((i-1)*(ep-1)+j,(i-1)*(ep-1)+k)+e_K(j,k)/e_j; %Multiplying by inverse jacobian in 1D comes down to just dividing by determinant
            end
        end
    end
    % Display the assembled mass and stiffness matrices
    M
    K
    
    %Apply clamped boundary condition at x=0 (remove first row and first column)
    M=M(2:end,2:end)
    K=K(2:end,2:end)
    
    % Compute eigenvalues and eigenvectors (natural frequencies and mode shapes)
    [eigenvect, eigenval] = eig(double(K),double(M))
    x_loc=linspace(0,l,p)
    hold off
    figure(800)
    plot(x_loc,-[0;eigenvect(:,1)]',x_loc,-[0;eigenvect(:,end)]')  % Plotting the first and last modal shape
    ylabel("First three modal shapes ")


    % Compute natural frequencies (square root of eigenvalues)
    omega_approx = sqrt(eigenval)
    
    % Convert frequency matrix to vector form
    omega_approx_vect=diag(omega_approx)
    first_natural_frequency = omega_approx(1)
    
    % Compute the exact natural frequencies for comparison
    omega=zeros(p-1,1);
    for k=1:p-1
        omega(k)=pi*(2*k-1)/2/l*sqrt(E/rho);
    end
    omega  % Display 
    first_exact_natural_frequency = omega(1)
    
    
    % Compute relative error in the predicted frequencies (factor of 4 error)
    relative_error=(omega_approx_vect-omega)./omega
    figure(1)
    plot(1:(p-1),relative_error)
    xlabel("Mode")
    ylabel("Relative error in predicted natural frequency")
    errfreq1=relative_error(1)  % Return the relative error in the first natural frequency
end