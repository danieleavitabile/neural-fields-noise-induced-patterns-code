clear all, close all, clc;

addpath('~/Dropbox/SecantContinuation/0.2/src');

%% Spatial coordinates
n = 50; 
x = linspace(0,1,n+2); y = x;
h = x(2) - x(1);
[X,Y] = meshgrid(x,x);
X = X(2:n+1,2:n+1); Y = Y(2:n+1,2:n+1);

%% Spatial operator
e = ones(n,1);
D2 = spdiags([e -2*e e], -1:1, n, n);
D2 = D2/h^2;

%% Laplacian
L = sparse( kron(D2,eye(n)) + kron(eye(n),D2) );

%% Initial parameters
p0(1) = 1;   % alpha
p0(2) = 0;   % lambda

%% Initial guess
u0 = zeros(n^2,1);

%% Assign problem 
problem = @(u,p) BratuProblem(u,p,L);

stepperPars.s0            = 0.5;
stepperPars.sMin          = 0.04;
stepperPars.sMax          = 3.0;
stepperPars.pMin          = 0;
stepperPars.pMax          = 8;
stepperPars.maxSteps      = 50;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 10;
stepperPars.iContPar      = 2;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'Jacobian','on');
stepperPars.optNonlinIter = 5;
stepperPars.dataFolder    = 'Data';
stepperPars.PlotSolution  = @(u,p,parentHandle) PlotSolution(u,p,parentHandle,X,Y,true);

branch = SecantContinuation(problem,u0,p0,stepperPars);
plot(branch(:,1),branch(:,2),'.-');
