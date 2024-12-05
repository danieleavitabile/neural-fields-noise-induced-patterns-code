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

%% Initial guess
% p0(1) = alpha;
% p0(2) = theta;
% p0(3) = D;
% p0(4) = A;
% u0 = cos(15*pi*x/L);

% sol = load('onebump.mat');
% u0  = sol.u;
% p0  = sol.p;

% sol = load('twobumps.mat');
% u0  = sol.u;
% p0  = sol.p;

% sol = load('threebumps.mat');
% u0  = sol.u;
% p0  = sol.p;

% sol = load('fourbumps.mat');
% u0  = sol.u;
% p0  = sol.p;

% sol = load('turing.mat');
% u0  = sol.u;
% p0  = sol.p;

% sol = load('turingunstable.mat');
% u0  = sol.u;
% p0  = sol.p;

% Homogeneous state
sol = load('HomogeneousSpectrum/solution_0000130.mat');
u0  = sol.u;
p0  = sol.p;

id = 2;
epsi = 10;
psi = real(sol.W(:,id));
u0 = u0 + epsi*psi;

% Plot initial solution
PlotSolution(x,u0,p0,[]);

% Define handle to right-hand side, jacobian and time output function
prob     = @(u,p) NeuralField(u,p,wHat,dx,[]);
jac      = @(u,p,v) NeuralFieldJacobianAction(u,p,wHat,dx,v);
plotSol  = @(u,p,parent) PlotSolution(x,u,p,parent);
solMeas  = @(step,u,p) SolutionMeasures(step,u,p,dx,2*L);
compSpec = [];%@(u,p) ComputeSpectrum(u,p,wHat,dx);
plotSpec = [];%@(lambda,p,parent) PlotSpectrum(lambda,p,parent);

%% Assign problem 
stepPars.iContPar                         = 3;
stepPars.pMin                             = -8;
stepPars.pMax                             = 8;
stepPars.s0                               =-0.003;
stepPars.sMin                             = 0.01;
stepPars.sMax                             = 0.5;
stepPars.maxSteps                         = 2000;
stepPars.nPrint                           = 1;
stepPars.nSaveSol                         = 1;
stepPars.finDiffEps                       = 1e-7;
stepPars.optNonlinIter                    = 5;
stepPars.NewtonGMRESOptions.nonlinTol     = 1e-5;
stepPars.NewtonGMRESOptions.nonlinMaxIter = 10;
stepPars.NewtonGMRESOptions.linTol        = 1e-2;
stepPars.NewtonGMRESOptions.linrestart    = 20;
stepPars.NewtonGMRESOptions.linmaxit      = 10;
stepPars.NewtonGMRESOptions.damping       = 1.0;
stepPars.NewtonGMRESOptions.display       = 0;
stepPars.dataFolder                       = 'Data';
stepPars.PlotSolution                     = plotSol;
stepPars.BranchVariables                  = solMeas;
stepPars.PlotBranchVariableId             = 1;
stepPars.ComputeEigenvalues               = compSpec;
stepPars.PlotSpectrum                     = plotSpec;

%% Run
branch = SecantContinuationNewtonGMRES(prob,jac,u0,p0,stepPars);

