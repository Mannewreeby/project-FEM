

T  = 20;
meshSizes = [1/50 1/20];

alpha1 = 0.01;
alpha2 = 0.005;
cf = 0.024;
ck = 0.055;
f = [cf; 0];
geometry = @circleg;







for i = 1:length(meshSizes)
   h = meshSizes(i);
   dt = h / 100;
   [p, e, t] = initmesh(geometry, 'hmax', h);
   [u, v] = initialConditionFunc(p, t);
   figure(1)
   pdesurf(p, t, u)
   figure(2)
   pdesurf(p, t, v)
   A = stiffnessAssembler2D(p, t);
   M = massAssembler2D(p, t);
   b = loadAssembler2D(p, t, f);
   
   
   
   for i = 1:T/dt
       S = u .* v.^2;
       u = (2 * M + dt*(cf*M - alpha1*A)) \ (2*dt*(b(:,1) - M*S) + (2 * M - dt*(cf*M - alpha1*A))*u);
       v = (2 * M + dt*((cf + ck)*M - alpha2*A)) \ (2*dt*(b(:, 2) + M*S) + (2 * M - dt*((cf + ck)*M - alpha2*A))*v);
       figure(1)
       pdesurf(p, t, u)
       title("u: " + num2str(i * dt))
       figure(2)
       pdesurf(p, t, v)
       title("v: " + num2str(i * dt))
       pause(0.1)
   end
   
   figure
   pdesurf(p, t, u);
   
   
   
   
end
