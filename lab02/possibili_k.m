clear
close all

kp=606;
kd=47;
g=9.8;
M=3;
m=0.1;
l=0.75
%X=Kd, Y=Kp
[X, Y] = meshgrid(0:500, [0:999]);

peaktime= Y >(1)/(4*M*l).* X.^(2)+(M+m)*g+pi^(2) * M*l; % Peak time condition
elong = Y < (pi^2+(log(5))^2)/(4*M*l*(log(5))^2)*X.^2+(M+m)*g; % Elongation condition
both = peaktime & elong; % Intersection of both inequalities

figure
colors = zeros(size(X)) + peaktime + elong; % Combine inequalities
scatter(X(:), Y(:), 3, colors(:), 'filled'); % Plot the region
hold on
scatter(kd,kp,'Marker','*');

