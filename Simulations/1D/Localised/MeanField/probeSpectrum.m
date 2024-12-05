%% Clean
clear all, close all, clc;

% addpath('~/Dropbox/SecantContinuation/0.2/src');
addpath('~/GitHub/neural-fields-interacting-particles/Simulations/SecantContinuation/0.2/src/');

%% Spatial grid
n = 2^10; L = 10*pi; dx = 2*L/n; x = -L+[0:n-1]'*dx;

%% Parameters
B = 0.4;

%% Funcion handles for the synaptic kernel
wFun = @(x) exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );
wHat = fft(wFun(x));

%% Solution 
% sol = load('fourbumps.mat');
% u0  = sol.u;
% p0  = sol.p;

% sol = load('TuringUnstable/solution_0000000.mat');
% u0  = sol.u;
% p0  = sol.p;

sol = load('Homogeneous/solution_0000026.mat');
u0  = sol.u;
p0  = sol.p;

% Plot initial solution
% PlotSolution(x,u0,p0,[]);

%% Funcion handles for the firing rate function
alpha = p0(1);
theta = p0(2);
D = p0(3);
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
S  = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));

W0 = integral(wFun, -inf,inf);
G = @(u) -u +W0*S(u);
uStar = 0.2;
uStar = fsolve(G,uStar)
u0 = uStar*ones(size(u0));

% Define handle to right-hand side, jacobian and time output function
prob     = @(u,p) NeuralField(u,p,wHat,dx,[]);
jac      = @(u,p,v) NeuralFieldJacobianAction(u,p,wHat,dx,v);
plotSol  = @(u,p,parent) PlotSolution(x,u,p,parent);
solMeas  = @(step,u,p) SolutionMeasures(step,u,p,dx,2*L);
compSpec = @(u,p) ComputeSpectrum(u,p,wHat,dx);
plotSpec = @(lambda,p,parent) PlotSpectrum(lambda,p,parent);

uStar = fsolve(@(u) prob(u,p0),u0);
[V,D] = compSpec(u0,p0);
lambda = diag(D);
plotSpec(lambda,p0,[]);

