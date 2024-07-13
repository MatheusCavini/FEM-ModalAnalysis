clear all; clc;

%Lista dos nós que compoem a estrutura
Nodes = {Node(0,0), Node(1.41,0), Node(0,4.1), Node(1.41,4.1), Node(0,6.4), ...
         Node(1.41,6.4), Node(0,7.25), Node(-1.5951, 8.135), ...
         Node(0.1122,9.1944), Node(1.6277,8.3194)};


%Função que calcula distância entre nós
function dist = Dist(nodeA, nodeB)
  dist = sqrt((nodeB.x - nodeA.x)^2 + (nodeB.y - nodeA.y)^2);
end

%Função que calcula o ângulo entre nós
function angle = Angle(nodeA, nodeB)
    angleRad = atan2(nodeB.y - nodeA.y, nodeB.x - nodeA.x);
    angle = rad2deg(angleRad);
end

%Propriedades do Material
E = 2.10*(10^11); %Pa
rho = 7600;     %Kg/m³

%Propriedades das Seções
A1 = 0.062;         %m²
A2 = 0.038;
A3 = 0.042;
A4 = 0.2199;
I4 = 0.013745;      %m^4
I3 = 0.0043814;

%Criação dos Pórticos da estrutura
portic1 = Portic(Dist(Nodes{1},Nodes{3}), E, A4, I4, rho, Angle(Nodes{1},Nodes{3}));
portic2 = Portic(Dist(Nodes{2},Nodes{4}), E, A4, I4, rho, Angle(Nodes{2},Nodes{4}));
portic3 = Portic(Dist(Nodes{3},Nodes{4}), E, A4, I4, rho, Angle(Nodes{3},Nodes{4}));
portic4 = Portic(Dist(Nodes{3},Nodes{5}), E, A3, I3, rho, Angle(Nodes{3},Nodes{5}));
portic5 = Portic(Dist(Nodes{4},Nodes{6}), E, A3, I3, rho, Angle(Nodes{4},Nodes{6}));

%Criação das Treliças da estrutura
link1 = Link(Dist(Nodes{5},Nodes{6}), E, A3, rho, Angle(Nodes{5},Nodes{6}));
link2 = Link(Dist(Nodes{5},Nodes{7}), E, A3, rho, Angle(Nodes{5},Nodes{7}));
link3 = Link(Dist(Nodes{6},Nodes{7}), E, A2, rho, Angle(Nodes{6},Nodes{7}));
link4 = Link(Dist(Nodes{7},Nodes{8}), E, A2, rho, Angle(Nodes{7},Nodes{8}));
link5 = Link(Dist(Nodes{8},Nodes{9}), E, A1, rho, Angle(Nodes{8},Nodes{9}));
link6 = Link(Dist(Nodes{7},Nodes{9}), E, A1, rho, Angle(Nodes{7},Nodes{9}));
link7 = Link(Dist(Nodes{9},Nodes{10}), E, A1, rho, Angle(Nodes{9},Nodes{10}));
link8 = Link(Dist(Nodes{6},Nodes{10}), E, A1, rho, Angle(Nodes{6},Nodes{10}));
link9 = Link(Dist(Nodes{7},Nodes{10}), E, A1, rho, Angle(Nodes{7},Nodes{10}));


% Listas dos elementos e conectividade de nós
elements = {portic1, portic2, portic3, portic4, portic5, link1, link2, link3, ...
            link4, link5, link6, link7, link8, link9};


connectivity = [1 3; 2 4; 3 4; 3 5; 4 6; 5 6; 5 7; 6 7; 7 8; 8 9; 7 9; 9 10; 6 10; 7 10];


% Montagem das matrizes globais de massa e rigidez
num_nodes = length(Nodes);
num_dof = 3 * num_nodes;   %3 DOFs por nó
K_global = zeros(num_dof);
M_global = zeros(num_dof);


for i = 1:length(elements)
    element = elements{i};
    nodes = connectivity(i, :);
    K_element = element.K;
    M_element = element.M;

    %Calcular os índices globais e explodir a matriz
     for r = 1:length(nodes)
        for s = 1:length(nodes)
            indices_r = (3*(nodes(r)-1)+1):(3*(nodes(r)-1)+3);
            indices_s = (3*(nodes(s)-1)+1):(3*(nodes(s)-1)+3);

            K_global(indices_r, indices_s) += ...
                K_element(3*(r-1)+1:3*(r-1)+3, 3*(s-1)+1:3*(s-1)+3);

            M_global(indices_r, indices_s) += ...
                M_element(3*(r-1)+1:3*(r-1)+3, 3*(s-1)+1:3*(s-1)+3);
        end
    end
end


%Aplicação das condições de contorno
function [K_constrained, M_constrained] = apply_constraints(K, M, constrained_dofs)
    % Get the total number of DOFs
    total_dofs = size(K, 1);

    % Create a logical mask for free DOFs
    free_dofs_mask = true(total_dofs, 1);
    free_dofs_mask(constrained_dofs) = false;

    % Apply the mask to the global matrices to get the constrained matrices
    K_constrained = K(free_dofs_mask, free_dofs_mask);
    M_constrained = M(free_dofs_mask, free_dofs_mask);
end

%Engasta os nós 1 e 2 e elimina rotação dos nós de treliça
constrained_dofs = [1, 2, 3, 4, 5, 6, 21, 24, 27, 30];

% Apply constraints
[K_constrained, M_constrained] = apply_constraints(K_global, M_global, constrained_dofs);



function [frequencies, mode_shapes] = modal_analysis(K_constrained, M_constrained)
    % Solve the generalized eigenvalue problem
    [V, D] = eig(K_constrained, M_constrained);

    % Extract eigenvalues and eigenvectors
    eigenvalues = diag(D);
    frequencies = sort(sqrt(eigenvalues) / (2 * pi));  % Convert to Hz
    mode_shapes = V;
end

[frequencies, mode_shapes] = modal_analysis(K_constrained, M_constrained);

% Display results
disp('Natural Frequencies (Hz):');
disp(frequencies(1:6));


