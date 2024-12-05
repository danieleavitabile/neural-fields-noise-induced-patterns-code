%% Cleaning
clear all, close all, clc;

%% Set path for input data
% dataPath = './';
dataPath = '/Users/daniele/GitHub/neural-fields-interacting-particles/Data/2D/Hexagon/Particles/';

%% Generate mesh and matrix data, or else load it; save animation
generateData = false;
saveAnimation = true;

%% Geometry and mesh parameters 
% R = 30; hmax = 0.01; 
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

% %% Model Parameters
% p(1) = 10;   % mu    
% p(2) = 0.6;  % theta 
% p(3) = 1;    % alpha
% p(4) = 2;    % beta
% p(5) = -2.2; % gamma
% p(6) = 1;    % delta
% p(7) = 3.5;  % nu     
% p(8) = 5;    % tau   
% p(9) = 0.2;  % D spiral
% % p(9) = 0.5;  % D travelling
% % p(9) = 1;  % D travelling
% % p(9) = 5;  % D travelling
% % p(9) = 10.0;  % D oscillations

%% Initial condition
x = mesh.nodes(:,1); y = mesh.nodes(:,2); n = length(x);
iU = 1:n; iV = n+iU;
% u0 = zeros(n,1); u0(find(y > 0)) = 1;
% v0 = zeros(n,1); v0(find(x < 0)) = 4;
% z0 = [u0; v0];
% u0 = zeros(n,1); u0(find(x > 0)) = 4;
% v0 = zeros(n,1); v0(find(x > 0)) = 4;
% z0 = [u0; v0];
sol  = load('spiral-h-0.03.mat');
z0 = sol.z; p = sol.p;
p(7) = 1.0;
% p(9) = 0.1;

%% Time step
D = p(9);
sigma = sqrt(2*D); ISparse = speye(2*n);
a = @(t,z) NeuralField(t,z,p,M);
b = @(t,z) sigma*ISparse;
tspan = [0 50];
dt = 0.1;
[~,~,t,zHist] = EulerMaruyama(a,b,tspan,z0,dt);

%% Animation
plotOpts.U.clim = [-2.5 2.5];
plotOpts.U.view = [0 90];
plotOpts.U.xlim = [-R R];
plotOpts.U.ylim = [-R R];
plotOpts.U.axis = 'equal';
plotOpts.U.markerSize = 10;
plotOpts.V = plotOpts.U;
PlotHistoryParticles(t,zHist,mesh.nodes,mesh.elements,plotOpts);

%% Save animation to file
if saveAnimation
  plotOpts.U.markerSize = 40;
  plotOpts.V.markerSize = 40;
  fileName = fullfile(dataPath,'animation.mp4');
  figure; SaveAnimationParticles(t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U,fileName);
end

figure;
plotOpts.U.markerSize = 30;
PlotAnimationSnapshot(1,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
exportgraphics(gcf,'t_0.0.eps','BackgroundColor','none','Resolution',400);
PlotAnimationSnapshot(51,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
exportgraphics(gcf,'t_5.0.eps','BackgroundColor','none','Resolution',400);
PlotAnimationSnapshot(101,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
exportgraphics(gcf,'t_10.0.eps','BackgroundColor','none','Resolution',400);
PlotAnimationSnapshot(500,t,zHist(:,iU),mesh.nodes,mesh.elements,plotOpts.U);
exportgraphics(gcf,'t_50.0.eps','BackgroundColor','none','Resolution',400);



