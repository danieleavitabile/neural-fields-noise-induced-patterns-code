clear all, close all, clc;

u = linspace(-5,5,1000);
mu = 20;
theta = 0.6;
S1 = @(u) 1./(1+exp(-mu*(u-theta)));

alpha = 50; D = 0.01;
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
S  = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));

plot(u,S1(u),u,S(u));
grid on;


