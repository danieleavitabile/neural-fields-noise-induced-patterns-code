%% Clean
clear all, close all, clc;

%% Spatial grid
n = 2^10; L = 10*pi; dx = 2*L/n; x = -L+[0:n-1]'*dx;

%% Parameters
% alpha = 10; theta  = 0.4; D = 10; B = 1.5; A = 30;
% alpha = 10; theta  = 0.4; D = 20; B = 1.5; A = 30;

% Increase noise to get patterns
% alpha = 10; theta  = 0.4; D = 0.02; B = 1.5; A = 7; % no patterns
alpha = 10; theta  = 0.4; D = 0.07; B = 1.5; A = 7;   % patterns

%% Funcion handles for the synaptic kernel
wFun = @(x) A/sqrt(pi)*exp(-x.^2)- A/(sqrt(pi)*B)*exp(-(x/B).^2); 
wHat = fft(wFun(x));

%% Funcion handles for the firing rate function
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
S  = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Plot Fourier Transform
xic = sqrt(8*log(B)/(B^2-1));

% %% Form the neural field RHS
p = [alpha; theta; D];
rhs = @(t,u) NeuralField(u,p,wHat,dx,[]);

%% Initial condition, and time step
% u0 = 0.3*rand(size(x));
u0 = 0.3*cos(xic*x);
tspan = [0 100];
[t,U] = ode45(rhs,tspan,u0);

%% Plot history
figure;
[X,T]= meshgrid(x,t);  
surf(X,T,U); view([0 90]); shading interp; caxis([-2 2]);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; colorbar;

p = [alpha; theta; D];
u = U(end,:)';
save('finalSate.mat','p','u');

