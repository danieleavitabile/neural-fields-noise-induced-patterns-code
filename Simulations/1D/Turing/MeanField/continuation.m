%% Clean
clear all, close all, clc;

% addpath('~/Dropbox/SecantContinuation/0.2/src');
addpath('~/GitHub/neural-fields-interacting-particles/Simulations/SecantContinuation/0.2/src/');

%% Spatial grid
n = 2^10; L = 10*pi; dx = 2*L/n; x = -L+[0:n-1]'*dx;

%% Parameters
% alpha = 10; theta  = 0.4; D = 10; B = 1.5; A = 30;
% alpha = 10; theta  = 0.4; D = 20; B = 1.5; A = 30;

% Increase noise to get patterns
% alpha = 10; theta  = 0.4; D = 0.02; B = 1.5; A = 7; % no patterns
alpha = 10; theta  = 0.4; D = 0.07; B = 1.5; A = 7;   % patterns

%% Funcion handles for the synaptic kernel
wFun = @(x) A/sqrt(pi)*exp(-x.^2)- A/(sqrt(pi)*B)*exp(-(x/B).^2); 
wHat = fft(wFun(x));

%% Funcion handles for the firing rate function
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
S  = @(u,D) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
% dS = @(u) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Initial guess
% p0(1) = alpha;
% p0(2) = theta;
% p0(3) = D;
% u0 = cos(15*pi*x/L);

sol = load('turing.mat');
u0  = sol.u;
p0  = sol.p;

u0 = 1*u0;
p0(3) = 0.001;

% Plot initial solution
% PlotSolution(x,u0,p0,[]);

% Define handle to right-hand side, jacobian and time output function
prob     = @(u,p) NeuralField(u,p,wHat,dx,[]);
jac      = @(u,p,v) NeuralFieldJacobianAction(u,p,wHat,dx,v);
plotSol  = @(u,p,parent) PlotSolution(x,u,p,parent);
solMeas  = @(step,u,p) SolutionMeasures(step,u,p,dx,2*L);
compSpec = @(u,p) ComputeSpectrum(u,p,wHat,dx);
plotSpec = @(lambda,p,parent) PlotSpectrum(lambda,p,parent);

%% Assign problem 
stepPars.iContPar                         = 3;
stepPars.pMin                             = 0;
stepPars.pMax                             = 8;
stepPars.s0                               = 0.003;
stepPars.sMin                             = 0.01;
stepPars.sMax                             = 0.1;
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

% %% Assign problem 
% stepPars.iContPar             = 3;
% stepPars.pMin                 = -8;
% stepPars.pMax                 = 8;
% stepPars.s0                   = 0.00003;
% stepPars.sMin                 =  0.0000001;
% stepPars.sMax                 =  0.00001;
% stepPars.maxSteps             = 20000;
% stepPars.nPrint               = 1;
% stepPars.nSaveSol             = 1;
% stepPars.finDiffEps           = 1e-7;
% stepPars.fsolveOptions        = optimset('Display','off',...
%                                          'TolFun',1e-4,...
%                                          'MaxIter',50,...
%                                          'Jacobian','off');
% stepPars.optNonlinIter        = 5;
% stepPars.dataFolder           = 'Data';
% stepPars.PlotSolution         = plotSol;
% stepPars.BranchVariables      = solMeas;
% stepPars.PlotBranchVariableId = 1;
% stepPars.ComputeEigenvalues   = compSpec;
% stepPars.PlotSpectrum         = plotSpec;

% %% Run
% branch = SecantContinuation(prob,u0,p0,stepPars);

