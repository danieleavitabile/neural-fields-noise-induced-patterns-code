%% Clean
clear all, close all, clc;

%% Parameters 
alpha = 10; theta  = 0.4; D = 0.07; B = 1.5; A = 7;   % patterns

%% Funcion handles for the synaptic kernel
wFun = @(x) A*(1/sqrt(pi)*exp(-x.^2)- 1/(sqrt(pi)*B)*exp(-(x/B).^2)); 

%% Firing rate
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
f = @(u) phi(alpha*(u-theta));

%% Number of particles and spatial grid
n = 2^10; L = 10*pi; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW

%%% CROSS CHECK SIGMA
sigma = sqrt(2*D);
wHat = fft(wFun(xi));
a = @(t,u) NeuralFieldFFT(u,wHat,f,dxi);
b = @(t,u) sigma*eye(length(u));

tspan = [0 50];
dt = 0.1;
u0 = zeros(n,1);
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

[XI,T] = meshgrid(xi,t);
surf(XI,T,U); shading interp; caxis([-2 2]); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; colorbar;

