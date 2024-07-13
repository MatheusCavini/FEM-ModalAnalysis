# Finite Elements Method: Modal Analysis
## Introduction
The code presented in this repository implements the functionalities needed to perform a modal analysis of mechanical structures made of trusses and beams in 2D space.<br>
The code was developed based on theorical concepts of Finite Elements Methods (FEM) presented in the discipline "Computational Mechanics" (PMR3401) as part of my Mechatronics Engineering course at Escola Politécnica da Universidade de São Paulo (USP).

## Classes
- **Node**: an object from the class `Node()` represents a node on 2D plane, that is part of a truss or beam. It takes only its X and Y coordinates.
- **Link**: an object from the class `Link()` represents a truss element. It takes its lenght `L`, Young's module `E`, area `A`, density `rho` and angle relative to the x axis `a`. On the backstages, it calculates the mass and stiffness matrices for the element.
- **Portic**: an object from the class `Portic()` represents a beam element. It takes its lenght `L`, Young's module `E`, area `A`, moment of inertia `I`, density `rho` and angle relative to the x axis `a`. On the backstages, it calculates the mass and stiffness matrices for the element.

## Usage example

Example on how to use the code provided to solve a modal analysis for the structure below, where 1 is a beam and 2 is a truss.
![image](https://github.com/user-attachments/assets/40ac5de0-133e-495f-ad95-e7a7efa4070a)

1. Instantiate the nodes that are part of the structure in a list:
```matlab
   %List of nodes
   Nodes = {Node(0,0), Node(2,0), Node(0, 2*sqrt(3))};
```

2. Define functions for calculating the distance and angle between 2 nodes:
```matlab
    %Function that calculates the distance between two nodes
    function dist = Dist(nodeA, nodeB)
      dist = sqrt((nodeB.x - nodeA.x)^2 + (nodeB.y - nodeA.y)^2);
    end
    
    %Function that calculates the angle between two nodes
    function angle = Angle(nodeA, nodeB)
      angleRad = atan2(nodeB.y - nodeA.y, nodeB.x - nodeA.x);
      angle = rad2deg(angleRad);
    end
```
3. Define the numerical properties for the elements (be consistent on units):
```matlab
  %Materials properties
  E = 10; %Pa
  rho1 = 3 ;
  rho2 = 2;     %Kg/m³
  
  %Sections properties
  A1 = 1;         %m²
  A2 = 2;
  I1 = 1;
```
4. Instantiate each element that is part of the structure:
```matlab
  %Instantiate elements
  portic1 = Portic(Dist(Nodes{1},Nodes{2}), E, A1, I1, rho1, Angle(Nodes{1},Nodes{2}));
  link1 = Link(Dist(Nodes{2},Nodes{3}), E, A2, rho2, Angle(Nodes{2},Nodes{3}));
```

5. Put the elements in a list. Create a connectivity matrix that obbey the order of the elements in the list, relating each one to its nodes:
```matlab
  %Elements list and connectivity matrix
  elements = {portic1, link1};
  connectivity = [1 2; 2 3];
```

6. Assemble the global mass and stiffness matrices from the local ones:
```matlab
  % Assemble the global mass and stiffness matrices
  num_nodes = length(Nodes);
  num_dof = 3 * num_nodes;   %3 DOFs por nó
  K_global = zeros(num_dof);
  M_global = zeros(num_dof);
  
  for i = 1:length(elements)
      element = elements{i};
      nodes = connectivity(i, :);
      K_element = element.K;
      M_element = element.M;
  
      % Compute global indexes and explode matrices
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
```

7. Apply the boundary conditions of the structure by removing the rows and columns of the matrices relative to the constrained degrees of freedom (DOF). Do the same to the 3rd degree of freedom of each truss-only node:
```matlab
  %Function to apply the boundary conditions
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
  
  %Encastre on node 1, support on 2 and revolute joint on 3 (besides, remove rotation of node 3 since is truss-only)
  constrained_dofs = [1, 2, 3, 5, 7, 8, 9];
  
  % Apply constraints
  [K_constrained, M_constrained] = apply_constraints(K_global, M_global, constrained_dofs);
```
8. Solve the Eigenvectors and Eigenvalues problem to the mass and stiffness matrices, and display the results:
```matlab
  function [frequencies, mode_shapes] = modal_analysis(K_constrained, M_constrained)
      % Solve the generalized eigenvalue problem
      [V, D] = eigs(K_constrained, M_constrained);
  
      % Extract eigenvalues and eigenvectors
      eigenvalues = diag(D);
      frequencies = sort(sqrt(eigenvalues));
      mode_shapes = V;
  end
  
  [frequencies, mode_shapes] = modal_analysis(K_constrained, M_constrained);
  
  % Display results
  disp('Natural Frequencies (ras/d):');
  disp(frequencies);
```



