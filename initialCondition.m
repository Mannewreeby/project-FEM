function [u, v] = initialCondition(p,t)
    % Function computing the initial condition according to the project
    % description
    
    % Number of elements
    nt = size(p,2);
    
    % Initiate initial condiation vector
    initial = zeros(length(p(1,:)),2);
    
    % Constants
    rho = 2;
    R = 0.5;
    r = 0.3;
    
    for K = 1:nt
        
        % Local-to-global mapping
        loc2glb = t(1:3,K);
        
        % Node x-coordinate
        x = p(1,loc2glb);
        
        % Node y-coordinate
        y = p(2,loc2glb);
        
        if abs(R-sqrt(x.^2+y.^2)) <= r
            initial(loc2glb, :) = [rho/2*rand() rho*(1 + 0.5*rand());rho/2*rand() rho*(1 + 0.5*rand());rho/2*rand() rho*(1 + 0.5*rand()) ];
        else
            initial(loc2glb, :) = [0 0; 0 0; 0 0];
        end
        
    end
    u = sparse(initial(:, 1));
    v = sparse(initial(:, 2));
    
end