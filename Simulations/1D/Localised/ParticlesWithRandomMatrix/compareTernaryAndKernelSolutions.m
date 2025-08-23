%% Clean
clear all, close all, clc;

%% Reproducible random simulation
rng(1,"twister");

%% Set true for using ternary matrix, false for matrix from kernel function
ternaryMatrix = true;

%% Set true for generating, flase for loading matrix
generateWeightMatrix = false;

%% Simulation type
% simType = 'Turing stable' 
% simType = 'Turing unstable'
simType = '3 bumps' 
switch simType 

  case '3 bumps'

    % Parameters
    alpha = 10; theta  = 0.9; D = 0.1; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.25; L = 10*pi; % 3 bumps

    % Function handle for initial condition
    u0Fun = @(x) u0Max./cosh(u0Alpha*x);

  case 'Turing stable'

    % Parameters
    alpha = 10; theta  = 0.9; D = 0.35; B = 0.4; A = 1.0; u0Max = 0.05; u0Alpha = 0.25; L = 10*pi; % Turing stable

    % Function handle for initial condition
    xic = 9*pi/L; u0Fun = @(x) u0Max*cos(xic*x);

  case 'Turing unstable'

    % Parameters
    alpha = 10; theta  = 0.9; D = 0.5; B = 0.4; A = 1.0; u0Max = 0.3; u0Alpha = 0.25; L = 10*pi; % Turing stable

    % Function handle for initial condition
    xic = 9*pi/L; u0Fun = @(x) u0Max*cos(xic*x);

  otherwise
    error('simType not known');

end

%% Funcion handle for the synaptic kernel
wFun = @(x) A*exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );

%% Funcion handles for the firing rates function (particles and mean field)
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
f = @(u) phi(alpha*(u-theta));

%% Number of particles and particles grid (their ratio must be a multiplee of 2)
n = 2^12; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Create sparse matrix
if ternaryMatrix
  if generateWeightMatrix 
    WMat = CreateWeightMatrix(wFun,xi);
    save('WMatrix.mat','WMat');
  else
    data =  load('WMatrix-n-4096.mat');
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

%% Function handle for the mean field
wHatXi = fft(wFun(xi));
rhs = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);

%% Initial condition (critical wavelength) and time span
dt = 0.1;
tspan = [0:dt:35];  % Here dt is used only for plotting, does not affect numerics
u0 = u0Fun(xi);

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
tspan = [0 70];
dt = 0.1; 
u0 = u0Fun(xi);
disp('Simulate Using Ternary Matrix')
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

%% Store ternary solutions
tTernary = t; UTernary = U;

%% Redo calculations with ternary
rng(1,"twister");
ternaryMatrix = false;

%% Function handle for the particle system
sigma = sqrt(2*D);
if ternaryMatrix % We use the ternary matrix
  a = @(t,u) -u + WMat*f(u);
else % We use the kernel
  wHatXi = fft(wFun(xi));
  a = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);
end
b = @(t,u) sigma*eye(length(u));

%% Function handle for the mean field
wHatXi = fft(wFun(xi));
rhs = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);

%% Initial condition (critical wavelength) and time span
dt = 0.1;
tspan = [0:dt:35];  % Here dt is used only for plotting, does not affect numerics
u0 = u0Fun(xi);

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
tspan = [0 70];
dt = 0.1; 
u0 = u0Fun(xi);
disp('Simulate Using Kernel')
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

%% Store solution with kernel
tKernel = t; UKernel = U;

%% Funcion handles for the firing rate function
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
S  = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Form the neural field RHS
p = [alpha; theta; D; A]; sigma = sqrt(2*D)
a = @(t,u) NeuralFieldMeanField(u,p,wHatXi,dxi,[]);
b = @(t,u) 0*eye(length(u));

u0 = u0Fun(xi);
dt = 0.1; 
tspan = [0 35];
disp('Simulate Mean Field')
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

%% Store mean field solution
tMeanField = t; UMeanField = U;

% pos = [2903 147 788 1117];
% cax = [-2 2.5];
% fh = figure;
% [XI,T] = meshgrid(xi,t);
% surf(XI,T,UTernary-UKernel); shading interp; caxis(cax); view([0 90]); 
% shading interp; xlabel('x'); ylabel('t'); zlabel('u^j_t');
% axis tight; % colorbar;
% axis off;
% set(fh,'position',pos);
% drawnow;
% exportgraphics(fh,'Figures/particle.jpeg','Resolution',300);

figure, hold on;
plot(xi,UTernary(end,:));
plot(xi,UKernel(end,:)); 
plot(xi, UMeanField(end,:), 'linewidth', 3);
xlim([-L L]); ylim([-2.5 4]); box;

