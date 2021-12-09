function M = massAssembler2D(p,t)
% Assembly of the Mass Matrix in 2D.

% Reference:
% Title: The Finite Element Method: Theory, Implementation, and Applications
% Author Mats G. Larson, Fredrik Bengzon
% Page: 63

%  Number of nodes
np = size(p,2);

% Number of elements
nt = size(t,2);

% Allocate mass matrix
M = zeros(np,np);

% Loop over all elements
for K = 1:nt
    % Local-to-global mapping
    loc2glb = t(1:3,K);
    % Node x-coordinate
    x = p(1,loc2glb);
    % Node y-coordinate
    y = p(2,loc2glb);
    % Area of triangle
    area = polyarea(x,y);
    % Element mass matrix
    MK = [2 1 1; 
          1 2 1;
          1 1 2;]/12*area;
    % Add element masses to m
    M(loc2glb,loc2glb) = M(loc2glb,loc2glb) + MK;
end
M = sparse(M);

end