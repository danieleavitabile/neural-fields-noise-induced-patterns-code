%% Clean
clear all, close all, clc;

%% Spatial grid
n = 2^10; L = 10*pi; dx = 2*L/n; x = -L+[0:n-1]'*dx;

%% Parameters
alpha = 10; theta  = 0.9; D = 0.5; B = 0.4; A = 1.0; u0Max = 0.3; u0Alpha = 0.25; L = 10*pi; % 3 bumps

%% Function handle for initial condition
xic = 9*pi/L; u0Fun = @(x) u0Max*cos(xic*x);

%% Funcion handles for the synaptic kernel
wFun = @(x) A*exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );
wHat = fft(wFun(x));

%% Funcion handles for the firing rate function
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
S  = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Form the neural field RHS
p = [alpha; theta; D; A]; sigma = sqrt(2*D)
a = @(t,u) NeuralField(u,p,wHat,dx,[]);
b = @(t,u) 0*eye(length(u));

%% Initial condition, and time step
% u0 = 0.3*rand(size(x));
% u0 = 0.3*cos(xic*x);
u0 = u0Fun(x);
dt = 0.1; 
tspan = [0 35];
 % [t,U] = ode45(a,tspan,u0);
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

%% Plot history
fh = figure;
pos = [2903 147 788 1117];
cax = [-1 2.2];
[X,T]= meshgrid(x,t);  
surf(X,T,U); view([0 90]); shading interp; caxis(cax);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; % colorbar;
axis off;
set(fh,'position',pos);
drawnow;
exportgraphics(fh,'meanfield.jpeg','Resolution',300);

p = [alpha; theta; D];
u = U(end,:)';
save('finalSate.mat','p','u');

