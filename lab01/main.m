clc; clear;

addpath(genpath('.'));

params.density = 1000;
params.h = 1;
params.d = 1;
params.g = 9.81;

pendulums = sloshing_pendulums(params, 10);
sys = sloshing_damped(pendulums, params, 0.003);

impulse(sys, 10);
grid on;
