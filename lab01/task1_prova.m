%% Laboratorio 1: Task 1
clc; clear;

density = 1000;
h = 1;
d = 1
g = 9.81;

M = d^2*pi/4*h*density;

n_modes = 10;

xi = zeros(n_modes,1)
xi(1:3) = [1.841; 5.329;  8.531]
for k = 4:n_modes
    xi(k) = xi(k-1) + pi
end

m = M * d * tanh(2*xi*h/d) ./ xi ./ (xi.^2-1) / h
m0 = M - sum(m)

L = d ./ (2 * xi .* tanh(2*xi*h/d))

H = h/2 - d/2./xi .* (tanh(xi*h/d) - 1 ./ sinh(2*xi*h/d))
H0 = sum(m.*(H-L))/m0

w_n = sqrt(g./L)
freq = w_n./(2*pi)

[["massa"; m0; m] ["lunghezza"; 0; L] ["distanza da CG"; H0; H] ["pulsazione"; 0; w_n] ["frequenza"; 0; freq]]