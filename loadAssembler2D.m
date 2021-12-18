function b = loadAssembler2D(p,t,f)
% Assembly of the Load vector.

% Reference:
% Title: The Finite Element Method: Theory, Implementation, and Applications
% Author Mats G. Larson, Fredrik Bengzon
% Page: 65

%  Number of nodes
np = size(p,2);

% Number of elements
nt = size(t,2);

% Allocate load matrix
b = zeros(np,1);

% Loop over elements
for K = 1:nt
    % Local-to-global mapping
    loc2glb = t(1:3,K);
    % Node x-cooridinate
    x = p(1,loc2glb);
    % Node y-cooridinate
    y = p(2,loc2glb);
    % Area
    area = polyarea(x,y);
    % Element load vector
    bK = [f(1);f(1);f(1)]/3*area;
    % Add element vector to b
    b(loc2glb) = b(loc2glb) + bK;  
end

end