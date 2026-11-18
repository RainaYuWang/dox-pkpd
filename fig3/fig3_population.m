function [base, density, A] = fig3_population(mesh, rate, kp, N0)
% With a=kp-rate, N=N0 exp(A)/(1+kp*N0*int(exp(A))/theta).
a = kp-rate;
increments = mesh.dt.*(a*mesh.weights);
A = [0; cumsum(increments)];
inner = A(1:end-1) + mesh.dt.*(a*mesh.partial);
B = [0; cumsum(mesh.dt.*(exp(inner)*mesh.weights))];
base = N0*exp(A(mesh.indices));
density = kp*N0*B(mesh.indices);
A = A(mesh.indices);
end
