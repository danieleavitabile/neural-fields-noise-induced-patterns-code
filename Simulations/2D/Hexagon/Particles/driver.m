%% Cleaning
clear all, close all, clc;

%% Set path for input data
% dataPath = './';
dataPath = '/Users/daniele/GitHub/neural-fields-interacting-particles/Data/2D/Hexagon/Particles/';

%% Generate mesh and matrix data, or else load it; save animation
generateData = false;
saveAnimation = false;

%% Geometry and mesh parameters 
% R = 30; hmax = 0.01; 
R = 30; hmax = 0.03; 

%% Model Parameters
p(1) = 10;   % mu    
p(2) = 0.6;  % theta 
p(3) = 1;    % alpha
p(4) = 2;    % beta
p(5) = -2.2; % gamma
p(6) = 1;    % delta
p(7) = 3.5;  % nu     
p(8) = 5;    % tau   
p(9) = 0.2;  % D spiral
% p(9) = 0.5;  % D travelling
% p(9) = 1;  % D travelling
% p(9) = 5;  % D travelling
% p(9) = 10.0;  % D oscillations

%% Load or generate mesh, synaptic matrix, and FEM matrix
if generateData

  % Mesh
  mesh = GenerateMesh(R,hmax,dataPath);

  % Synaptic matrix 
  r = linspace(-60,60,1000);
  w = SynapticKernelBessel(r);
  wFun = griddedInterpolant(r,w);
  W = GenerateSynapticMatrix(mesh.nodes,wFun,1e-3,dataPath);

  % FEM matrix
  M = GenerateFEMMatrix(mesh,W,dataPath);

else

  % Mesh
  fileName = fullfile(dataPath,'mesh-h-0.03.mat');
  mesh = load(fileName);

  % Synaptic matrix
  fileName = fullfile(dataPath,'synaptic-matrix-h-0.03.mat');
  data = load(fileName); W = data.W;

  % FEM matrix
  fileName = fullfile(dataPath,'fem-matrix-h-0.03.mat');
  data = load(fileName); M = data.M;

end

%% Initial condition
x = mesh.nodes(:,1); y = mesh.nodes(:,2); n = length(x);
iU = 1:n; iV = n+iU;
u0 = zeros(n,1); u0(find(y > 0)) = 1;
v0 = zeros(n,1); v0(find(x < 0)) = 4;
z0 = [u0; v0];
% u0 = zeros(n,1); u0(find(x > 0)) = 4;
% v0 = zeros(n,1); v0(find(x > 0)) = 4;
% z0 = [u0; v0];
% sol  = load('spiral-h-0.03.mat');
% z0 = sol.z; p = sol.p;
% % p(7) = 1.0;

%% Time step
D = p(9);
sigma = sqrt(2*D); ISparse = speye(2*n);
a = @(t,z) NeuralField(t,z,p,M);
b = @(t,z) sigma*ISparse;
tspan = [0 50];
dt = 0.1;
[~,~,t,zHist] = EulerMaruyama(a,b,tspan,z0,dt);
% [t,zHist] = ode45(rhs,tspan,z0);

%% Animation
plotOpts.U.clim = [-2 5];
plotOpts.U.view = [0 90];
plotOpts.U.xlim = [-R R];
plotOpts.U.ylim = [-R R];
plotOpts.U.axis = 'equal';
plotOpts.U.markerSize = 10;
plotOpts.V = plotOpts.U;
PlotHistoryParticles(t,zHist,mesh.nodes,mesh.elements,plotOpts);

%% Save animation to file
if saveAnimation
  fileName = fullfile(dataPath,'animation.mp4');
  figure; SaveAnimationParticles(t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U,fileName);
end

%% Reset initial condition and restart 
u0 = zHist(end,iU)' .* ( sqrt(x.^2 + y.^2) > 15 );
v0 = zHist(end,iV)' .* ( sqrt(x.^2 + y.^2) > 15 );
z0 = [u0; v0];
tspan = [0 50];
dt = 0.1;
[t,zHist] = EulerMaruyama(a,b,tspan,z0,dt);

% Animation
PlotHistoryParticles(t,zHist,mesh.nodes,mesh.elements,plotOpts);

% Save animation to file
if saveAnimation
  fileName = fullfile(dataPath,'animation-perturbed.mp4');
  figure; SaveAnimationParticles(t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U,fileName);
end
