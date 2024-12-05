clear all, close all, clc;

p0=[1; 0];
u0=[1; 0];

stepperPars.s0            = -0.1;
stepperPars.maxSteps      = 30;
stepperPars.iContPar      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','none',...
                                     'Jacobian','on');

problem = @(u,p) FoldProblem(u,p);
SecantContinuation(problem,u0,p0,stepperPars);
