function integral = integralApproximation(p,t,xi)
% Approximate integral in 2D using the Trapezodial rule

% Number of elements
nt = size(t,2);

% Initialize the integral
integral = 0;

for K = 1:nt
    
    % Local-to-global mapping
    loc2glb = t(1:3,K);
    % Node x-coordinate
    x = p(1,loc2glb);
    
    % Node y-coordinate
    y = p(2,loc2glb);
    
    % Calculate area
    area = polyarea(x,y);
    % Summation
    sum = 0;
    for j = 1:length(loc2glb)
        sum = sum + xi(loc2glb(j));
    end
    % Add result to integral
    integral = integral + sum*area/3;
end

end