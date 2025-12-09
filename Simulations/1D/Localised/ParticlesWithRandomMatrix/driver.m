%% Clean
clear all, close all, clc;

%% Reproducible random simulation
rng(1,"twister");

%% Set true for using ternary matrix, false for matrix from kernel function
ternaryMatrix = false;

%% Set true for generating, flase for loading matrix
generateWeightMatrix = false;

%% Simulation type
% simType = 'Turing stable' 
simType = 'Turing unstable'
% simType = '3 bumps' 
switch simType 

  case '3 bumps'

    % Parameters
    alpha = 10; theta  = 0.9; D = 0.1; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.25; L = 10*pi; 

    % Function handle for initial condition
    u0Fun = @(x) u0Max./cosh(u0Alpha*x);

  case 'Turing stable'

    % Parameters
    alpha = 10; theta  = 0.9; D = 0.35; B = 0.4; A = 1.0; u0Max = 0.3; u0Alpha = 0.25; L = 10*pi;

    % Function handle for initial condition
    xic = 9*pi/L; u0Fun = @(x) u0Max*cos(xic*x);

  case 'Turing unstable'

    % Parameters
    alpha = 10; theta = 0.9; D = 0.5; B = 0.4; A = 1.0; u0Max = 0.3; u0Alpha = 0.25; L = 10*pi; 

    % Function handle for initial condition
    xic = 9*pi/L; u0Fun = @(x) u0Max*cos(xic*x);

  otherwise
    error('simType not known');

end
sigma = sqrt(2*D)

%% Funcion handle for the synaptic kernel
wFun = @(x) A*exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );


%% Funcion handles for the firing rates function (particles and mean field)
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
f = @(u) phi(alpha*(u-theta));

%% Number of particles and particles grid (their ratio must be a multiplee of 2)
n = 2^12; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

% figure, hold on;
% plot(xi,  0.05 + wFun(xi).* ( wFun(xi) > 0))
% plot(xi, -0.05 + wFun(xi).* ( wFun(xi) < 0))
% % plot(xi, wFun(xi))
% axis tight; 
% axis off; 
% pause

%% Create sparse matrix
if ternaryMatrix
  if generateWeightMatrix 
    WMat = CreateWeightMatrix(wFun,xi);
    save('WMatrix.mat','WMat');
  else
    data =  load('WMatrix-n-128.mat');
    % data =  load('WMatrix-n-4096.mat');
    % data =  load('WMatrix-n-16384.mat');
    WMat = data.WMat;
  end
end

%% Function handle for the particle system
sigma = sqrt(2*D);
if ternaryMatrix % We use the ternary matrix
  a = @(t,u) -u + WMat*f(u);
else % We use the kernel
  wHatXi = fft(wFun(xi));
  a = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);
end
b = @(t,u) sigma*eye(length(u));

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
tspan = [0 35];
dt = 0.1; 
u0 = u0Fun(xi);
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

pos = [2903 147 788 1117];
cax = [-1 2.2];
fh = figure;
[XI,T] = meshgrid(xi,t);
surf(XI,T,U); shading interp; caxis(cax); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u^j_t');
axis tight; % colorbar;
axis off;
set(fh,'position',pos);
drawnow;
exportgraphics(fh,'Figures/particle.jpeg','Resolution',300);

