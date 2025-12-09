%% Clean
clear all, close all, clc;

%% Parameters 
% alpha = 10; theta  = 0.9; D = 0.0; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.25; L = 10*pi; % 1 bump
alpha = 10; theta  = 0.9; D = 0.1; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.25; L = 10*pi; % 3 bumps
% alpha = 10; theta  = 0.9; D = 0.5; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.25; L = 10*pi; % Turing

%% Funcion handle for the synaptic kernel
wFun = @(x) A*exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );

%% Funcion handles for the firing rates function (particles and mean field)
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
f = @(u) phi(alpha*(u-theta));

%% Function handle for initial condition
u0Fun = @(x) u0Max./cosh(u0Alpha*x);

%% Number of particles and particles grid (their ratio must be a multiplee of 2)
n = 2^12; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Kerenel's FFT in both mean field and particle system
wHatXi = fft(wFun(xi));

%% Function handle for the particle system
sigma = sqrt(2*D);
a = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);
b = @(t,u) sigma*eye(length(u));

%% Function handle for the mean field
rhs = @(t,u) NeuralFieldFFT(u,wHatX,S,dx);

%% Initial condition (critical wavelength) and time span
dt = 0.1;
tspan = [0:dt:35];  % Here dt is used only for plotting, does not affect numerics
u0 = u0Fun(xi);

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
tspan = [0 35];
dt = 0.1; 
u0 = u0Fun(xi);
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

pos = [2903 147 788 1117];
cax = [-2 2.5];
fh = figure;
[XI,T] = meshgrid(xi,t);
surf(XI,T,U); shading interp; caxis(cax); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u^j_t');
axis tight; % colorbar;
axis off;
set(fh,'position',pos);
drawnow;
exportgraphics(fh,'Figures/particle.jpeg','Resolution',300);

