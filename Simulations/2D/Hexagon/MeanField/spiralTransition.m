%% Cleaning
clear all, close all, clc;

%% Set path for input data
% dataPath = './';
dataPath = '/Users/daniele/GitHub/neural-fields-noise-induced-patterns-code/Data/2D/Hexagon/MeanField/';

%% Generate mesh and matrix data, or else load it; save animation
generateData = false;
saveAnimation = true;

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
% The .mat file has D = sigma = 0, uncomment here for D = 0.1, that is, sigma = 0.4472 
p(9) = 0.1; 

%% Time step
rhs = @(t,u) NeuralField(t,u,p,M);
tspan = [0 100];
% [t,zHist] = ode45(rhs,tspan,z0);
b = @(t,u) 0*t;
dt = 0.1;
[~,~,t,zHist] = EulerMaruyama(rhs,b,tspan,z0,dt);

%% Animation
plotOpts.U.clim = [-2.5 2.5];
plotOpts.U.view = [0 90];
plotOpts.U.xlim = [-R R];
plotOpts.U.ylim = [-R R];
plotOpts.U.axis = 'equal';
plotOpts.U.markerSize = 10;
plotOpts.V = plotOpts.U;
PlotHistory(t,zHist,mesh.nodes,mesh.elements,plotOpts);

%% Save animation to file
if saveAnimation
  plotOpts.U.markerSize = 40;
  plotOpts.V.markerSize = 40;
  fileName = fullfile(dataPath,'animation.mp4');
  figure; SaveAnimation(t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U,fileName);
end

% figure;
% plotOpts.U.markerSize = 30;
% PlotAnimationSnapshot(1,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
% exportgraphics(gcf,'t_0.0.eps','BackgroundColor','none','Resolution',400);
% PlotAnimationSnapshot(51,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
% exportgraphics(gcf,'t_5.0.eps','BackgroundColor','none','Resolution',400);
% PlotAnimationSnapshot(101,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
% exportgraphics(gcf,'t_10.0.eps','BackgroundColor','none','Resolution',400);
% PlotAnimationSnapshot(500,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
% exportgraphics(gcf,'t_50.0.eps','BackgroundColor','none','Resolution',400);
figure;
plotOpts.U.markerSize = 30;
PlotAnimationSnapshot(1,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
exportgraphics(gcf,'t_0.0.eps','BackgroundColor','none','Resolution',400);
PlotAnimationSnapshot(11,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
exportgraphics(gcf,'t_1.0.eps','BackgroundColor','none','Resolution',400);
PlotAnimationSnapshot(51,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
exportgraphics(gcf,'t_5.0.eps','BackgroundColor','none','Resolution',400);
PlotAnimationSnapshot(501,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
exportgraphics(gcf,'t_50.0.eps','BackgroundColor','none','Resolution',400);

fileName = fullfile(dataPath,'historyData.mat');
save(fileName, 't', 'zHist', 'mesh');

