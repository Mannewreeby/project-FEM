function A = stiffnessAssembler2D(p,t)
% Assembly of the Stiffness Matrix in 2D

% Reference:
% Title: The Finite Element Method: Theory, Implementation, and Applications
% Author Mats G. Larson, Fredrik Bengzon
% Page: 86

%  Number of nodes
np = size(p,2);

% Number of elements
nt = size(t,2);

% Allocate stiffness matrix
A = zeros(np,np);

% Loop over elements
for K = 1:nt
    % Local to global map
    loc2glb = t(1:3,K); 
    % Node x-coordinate
    x = p(1,loc2glb); 
    % Node y-coordinate
    y = p(2,loc2glb); 
    
    % Fetch the gradients of the hat functions
    [area,b,c] = hatGradients(x,y);
    
    % Element stiffness matrix
    
    AK = (b*b.' + c*c.')*area; 
    
    %Add element stiffnesses to A
    A(loc2glb, loc2glb) = A(loc2glb, loc2glb) + AK;
end
A = sparse(A);

end
