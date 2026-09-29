clear
close all
clc

%engenering wood
% E: modulus of elasticity
E = 12e9;
% A: area of cross section
A = 100e-4; 
% rho: density
rho = 600; 

% rhoA: mass per unit length 
rhoA = rho*A;
% EA: axial stiffness
EA = E*A; 

multi = 10; % multiplication factor for displacement plot


syms x 


% Main variable 1
nb_node_inter=4;


% construction of the node and connectivity variable
connectivity=[1 2;1 6;2 6;6 7;1 7];
for inter_node_i=3:2:(nb_node_inter-1)
    real_node_index=inter_node_i+5;
    connectivity=[connectivity;real_node_index real_node_index-2;real_node_index real_node_index-1;real_node_index real_node_index+1;real_node_index-1 real_node_index+1];
end
real_node_index=nb_node_inter+5;
connectivity=[connectivity;real_node_index 5;real_node_index-1 5;real_node_index real_node_index+1;real_node_index+1 5;
real_node_index+2 5;real_node_index+1 real_node_index+2];
for inter_node_i=nb_node_inter+3:2:(2*nb_node_inter-1)
    real_node_index=inter_node_i+5;
    connectivity=[connectivity;real_node_index real_node_index-2;real_node_index real_node_index-1;real_node_index-1 real_node_index+1;real_node_index real_node_index+1];
end
real_node_index=2*nb_node_inter+5;
connectivity=[connectivity;real_node_index 4;real_node_index 3;real_node_index-1 3;3 4];
 
%main variable 2
a=2.5;

%shape function
f(x) = -a/121*(x+1)*(x-21);
g(x)= f(x)-f(0);


nodesCoordinates = [ 0 0
                    -1 0
                    20 0
                    21 0
                    10 a];
node_spacing=10/nb_node_inter;
for inter_node_i=1:(nb_node_inter*2)
    if inter_node_i<=nb_node_inter
        x_node=(inter_node_i-1)*node_spacing;
        if rem(inter_node_i, 2) == 0
            nodesCoordinates=[nodesCoordinates;x_node g(x_node)];
        else
            nodesCoordinates=[nodesCoordinates;x_node f(x_node)];
        end
    else
        x_node=(inter_node_i-1)*node_spacing+node_spacing;
        if rem(inter_node_i, 2) == 0
            nodesCoordinates=[nodesCoordinates;x_node f(x_node)];
        else
            nodesCoordinates=[nodesCoordinates;x_node g(x_node)];
        end
    end
end
    

%solve the bridge
numberOfNodes = size(nodesCoordinates,1);
GDof = 2*numberOfNodes;

% compute and assemble the structure stiffness matrix
stiffness = formStiffness2Dtruss(GDof, connectivity, nodesCoordinates, EA);

% compute and assemble the structure maass matrix
consistentMass = formConsistentMass2Dtruss(GDof,connectivity, nodesCoordinates, rhoA);
lumpedMass = formLumpedMass2Dtruss(GDof,connectivity, nodesCoordinates, rhoA);

prescribedDof = transpose([1 2 3 4 5 6 7 8]);


[modes_consisMass,omega_consisMass] = computeFrequenciesAndModes(GDof,prescribedDof,stiffness,consistentMass,0);
[modes_lumpedMass,omega_lumpedMass] = computeFrequenciesAndModes(GDof,prescribedDof,stiffness,lumpedMass,0);

% Extract the first 3 natural frequencies (in Hz)
freq_consis = omega_consisMass(1:3)/(2*pi);
freq_lumped = omega_lumpedMass(1:3)/(2*pi);
FrequencyTable = table((1:3)', freq_consis, freq_lumped, ...
    'VariableNames', {'ModeNumber', 'ConsistentMass_Hz', 'LumpedMass_Hz'});
display(FrequencyTable)

%remove the displacement of 0 of the stiffness matrix
stiffness_reduced=stiffness;
stiffness_reduced(prescribedDof,:)=[];
stiffness_reduced(:,prescribedDof)=[];

GDof_reduced=GDof-length(prescribedDof);


%get the dispalcement F=Kx -> x=F/K
F=zeros(GDof_reduced,1);
F(2)=-50000;

%displacement dx_tot as a line and dx_nodes in the same shape as
%nodecoordinates
dx=stiffness_reduced\F;
dx_tot=zeros(1,2*numberOfNodes);
dx_tot(prescribedDof)=zeros(1,8);
dx_tot(9:2*numberOfNodes)=dx';
dx_nodes=reshape(dx_tot,2,numberOfNodes)';

% get the lengths of the bar

numBars = size(connectivity, 1);% Number of bars

% Initialize length vectors
L0 = zeros(numBars, 1);
Ld = zeros(numBars, 1);
stess=zeros(numBars, 1);
strain=zeros(numBars, 1);

% Compute deformed positions
nodes_deformed = nodesCoordinates + dx_nodes;

% Compute exagerated deformed positions with the multi parameter
nodes_deformed_exagerated = nodesCoordinates + multi*dx_nodes;

for i = 1:numBars
    % Get the indices of the nodes forming the bar
    n1 = connectivity(i, 1);
    n2 = connectivity(i, 2);
    
    % Compute initial length
    L0(i) = norm(nodesCoordinates(n2, :) - nodesCoordinates(n1, :));
    
    % Compute deformed length
    Ld(i) = norm(nodes_deformed(n2, :) - nodes_deformed(n1, :));

    % Compute strain
    strain(i) = (Ld(i) - L0(i)) / L0(i);
    
    % Compute stress using Hooke's Law
    stress(i) = E * strain(i);
end
delta_L_table=table((1:numBars)',(Ld-L0)*1000,strain*100,stress'/10^6,'VariableNames',{'Bar Number','Deformation in mm','strain [%]','stress [MPa]'});
display(delta_L_table)
display('volume of the beam:'+string(sum(L0*A))+'m^2')




% % Set LaTeX as the interpreter for text in plots
% set(groot, 'defaultTextInterpreter', 'latex');
% set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% set(groot, 'defaultLegendInterpreter', 'latex');

% to check if the connectivity and nodecoordinates are correct
% figure
% draw2Dtruss(nodesCoordinates, connectivity, ...
%     'LineColor', 'r', ...
%     'LineWidth', 2, ...
%     'MarkerSize', 2, ...
%     'ShowNodeNumbers', true,'ShowBeamNumbers',true);


% mode number visu
for modeNumber=[1,2,3]
    
    us = 1:2:2*numberOfNodes-1;
    vs = 2:2:2*numberOfNodes;
    
    full_modes = zeros(GDof,length(omega_consisMass));
    activeDof = setdiff(transpose((1:GDof)), prescribedDof);
    
    full_modes(activeDof,:) = modes_consisMass;
    
    XX = full_modes(us,modeNumber); 
    YY = full_modes(vs,modeNumber);
    dispNorm = max(sqrt(XX.^2+YY.^2));
    scaleFact = 2*dispNorm;
    
    figure('Position',[100 100 600 200])
    draw2Dtruss(nodesCoordinates, connectivity, ...
        'LineColor', [1 0 0 .2], ...
        'LineWidth', 2, ...
        'MarkerSize', 2, ...
        'ShowNodeNumbers', false,'ShowBeamNumbers',false);
    draw2Dtruss(nodesCoordinates+scaleFact*[XX YY], connectivity, ...
        'LineColor', 'r', ...
        'LineWidth', 2, ...
        'MarkerSize', 2, ...
        'ShowNodeNumbers', true,'ShowBeamNumbers',false,'Title','Mode Shape '+string(modeNumber));
    exportgraphics(gcf,'mode'+string(modeNumber)+'.pdf',"ContentType","vector")
end

%plot the example figure
figure('Position',[100 100 600 200])
draw2Dtruss(nodesCoordinates, connectivity, ...
    'LineColor', [1 0 0], ...
    'LineWidth', 2, ...
    'MarkerSize', 15, ...
    'ShowNodeNumbers', true,'ShowBeamNumbers',false);
x_up=linspace(-1,21,500);
plot(x_up,f(x_up),'LineWidth',1,'LineStyle','--','color',[0 0 0])
hold on
x_down=linspace(0,20,500);
plot(x_down,g(x_down),'LineWidth',1,'LineStyle','--','color',[0 0 0])
exportgraphics(gcf,'shapefig.pdf','ContentType','vector')

%plot the table of the displacement of the nodes
table_displacement=table((5:13)',string(round(sqrt(dx_nodes(5:13,1).^2+dx_nodes(5:13,2).^2)*1000,2))+' mm',VariableNames={'nodes','displacement'});
display(table_displacement);

%plot the figure with the axial load colormap
colormapplot(nodesCoordinates,connectivity,dx_nodes,nodes_deformed_exagerated,E)
exportgraphics(gcf,'colormap_opti.pdf','ContentType','vector')
