%% Cleaning
clear all, close all, clc;

%% Set path for input data
% dataPath = './';
dataPath = '/Users/daniele/GitHub/neural-fields-interacting-particles/Data/2D/Hexagon/MeanField/';

%% Generate mesh and matrix data, or else load it; save animation
generateData = false;
saveAnimation = false;

%% Geometry and mesh parameters 
R = 30; hmax = 0.03; 

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
sol  = load('spiral-h-0.03.mat');
z0 = sol.z; p = sol.p;
p(7) = 1.0;
% p(9) = 0.1;

%% Time step
rhs = @(t,u) NeuralField(t,u,p,M);
tspan = [0 50];
[t,zHist] = ode45(rhs,tspan,z0);

%% Animation
plotOpts.U.clim = [-2.5 2.5];
plotOpts.U.view = [0 90];
plotOpts.U.xlim = [-R R];
plotOpts.U.ylim = [-R R];
plotOpts.U.axis = 'equal';
plotOpts.V = plotOpts.U;
PlotHistory(t,zHist,mesh.nodes,mesh.elements,plotOpts);

%% Save animation to file
if saveAnimation
  fileName = fullfile(dataPath,'animation.mp4');
  figure; SaveAnimation(t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U,fileName);
end

