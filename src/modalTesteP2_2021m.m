clear all; clc;

%Lista dos nós que compoem a estrutura
Nodes = {Node(0,0), Node(2,0), Node(0, 3.4641)};


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
E = 10; %Pa
rho1 = 3 ;
rho2 = 2;     %Kg/m³

%Propriedades das Seções
A1 = 1;         %m²
A2 = 2;
I1 = 1;


%Criação dos Pórticos da estrutura
portic1 = Portic(Dist(Nodes{1},Nodes{2}), E, A1, I1, rho1, Angle(Nodes{1},Nodes{2}));
link1 = Link(Dist(Nodes{2},Nodes{3}), E, A2, rho2, Angle(Nodes{2},Nodes{3}));



% Listas dos elementos e conectividade de nós
elements = {portic1, link1};

connectivity = [1 2; 2 3];

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

%Engaste em 1, articulação em 3 e apoio em 2
constrained_dofs = [1, 2, 3, 5, 7, 8, 9];

% Apply constraints
[K_constrained, M_constrained] = apply_constraints(K_global, M_global, constrained_dofs);

function [frequencies, mode_shapes] = modal_analysis(K_constrained, M_constrained)
    % Solve the generalized eigenvalue problem
    [V, D] = eigs(K_constrained, M_constrained);

    % Extract eigenvalues and eigenvectors
    eigenvalues = diag(D);
    frequencies = sort(sqrt(eigenvalues) );  % Convert to Hz
    mode_shapes = V;
end

[frequencies, mode_shapes] = modal_analysis(K_constrained, M_constrained);

% Display results
disp('Natural Frequencies (ras/d):');
disp(frequencies);
