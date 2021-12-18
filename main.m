close all;
T = 20;
meshSizes = [1/5 1/20];

alpha1 = 0.01;
alpha2 = 0.005;
cf = 0.024;
ck = 0.055;
f = [cf; 0];
geometry = @circleg;

for i = 1:length(meshSizes)
    h = meshSizes(i);
    dt = h / 2;
    [p, e, t] = initmesh(geometry, 'hmax', h);
    [xi1Initial, xi2Initial] = initialCondition(p, t);
    xi1 = xi1Initial; xi2 = xi2Initial;


    figure()
    title('u initial, h=' + string(h))
    pdesurf(p, t, xi1)
    xlabel("x")
    ylabel('y')
    zlabel('u')


    figure()
    title('v initial, h=' + string(h))
    pdesurf(p, t, xi2)
    xlabel("x")
    ylabel('y')
    zlabel('v')


    A = stiffnessAssembler2D(p, t);
    M = massAssembler2D(p, t);
    b = loadAssembler2D(p, t, f);

    massLoss_u = zeros(T / dt, 1);
    massLoss_v = zeros(T / dt, 1);

    time = 0:dt:T;
    k = 1;

    for j = 1:length(time)

        if j > 1
            S = xi1 .* xi2.^2;
            xi1 = (2 * M + dt * (cf * M + alpha1 * A)) \ (2 * dt * (b - M * S) + (2 * M - dt * (cf * M + alpha1 * A)) * xi1);
            xi2 = (2 * M + dt * ((cf + ck) * M + alpha2 * A)) \ (2 * dt * (M * S) + (2 * M - dt * ((cf + ck) * M + alpha2 * A)) * xi2);
        end

        %        figure(1)
        %        pdesurf(p, t, xi1)
        %        title("u: " + num2str(time))
        %        figure(2)
        %        pdesurf(p, t, xi2)
        %        title("v: " + num2str(time))
        g_u = xi1Initial - xi1;
        g_v = xi2Initial - xi2;
        massLoss_u(j) = integralApproximation(p, t, g_u);
        massLoss_v(j) = integralApproximation(p, t, g_v);

    end

    figure()
    pdesurf(p, t, xi1);
    title('u final, h=' + string(h))
    xlabel("x")
    ylabel('y')
    zlabel('u')

    figure()
    pdesurf(p, t, xi2);
    title('v final, h=' + string(h))
    xlabel("x")
    ylabel('y')
    zlabel('v')

    figure()
    plot(time, massLoss_u);
    title('Mass loss u, h=' + string(h))
    xlabel('time')
    ylabel('mass loss')

    figure()
    plot(time, massLoss_v);
    title('Mass loss v, h=' + string(h))
    xlabel('time')
    ylabel('mass loss')

end
