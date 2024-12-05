%% Clean
clear all, close all, clc;

%% Parameters 
alpha = 10; theta  = 0.4; D = 0.17; B = 1.5; A = 7; L = 10*pi;% patterns

%% Funcion handle for the synaptic kernel
wFun = @(x) A*(1/sqrt(pi)*exp(-x.^2)- 1/(sqrt(pi)*B)*exp(-(x/B).^2)); 

%% Funcion handles for the firing rates function (particles and mean field)
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
f = @(u) phi(alpha*(u-theta));
S = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));

%% Function handle for initial condition
lambdac = sqrt(8*log(B)/(B^2-1)); 
u0Fun = @(x) 0.01*cos(lambdac*x);

%% Mean field grid
q = 2^7; dx = 2*L/q; x = -L+[0:q-1]'*dx;

%% Number of particles and particles grid (their ratio must be a multiplee of 2)
r = 2^0; n = q*r; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Kerenel's FFT in both mean field and particle system
wHatXi = fft(wFun(xi));
wHatX  = fft(wFun(x));

%% Function handle for the particle system
sigma = sqrt(2*D);
a = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);
b = @(t,u) sigma*eye(length(u));

%% Function handle for the mean field
rhs = @(t,u) NeuralFieldFFT(u,wHatX,S,dx);

%% Initial condition (critical wavelength) and time span
dt = 0.1;
tspan = [0:dt:35];  % Here dt is used only for plotting, does not affect numerics
u0 = u0Fun(x);

%% Timestep mean field
[tmf,Umf] = ode45(rhs,tspan,u0);

%% Plot
pos = [2903 147 788 1117];
cax = [-2.5 2.5];
figure;
% subplot(2,2,[1,3]);
[X,T] = meshgrid(x,tmf);
surf(X,T,Umf); caxis(cax); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; colorbar;
axis off;
set(gcf,'position',pos);
drawnow;
exportgraphics(gcf,'Figures/meanfield.jpeg','Resolution',300);
pause

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
tspan = [0 35];
dt = 0.1; 
u0 = u0Fun(xi);
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

% subplot(2,2,[2,4]);
[XI,T] = meshgrid(xi,t);
surf(XI,T,U); shading interp; caxis(cax); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; % colorbar;
set(gcf,'position',pos);
drawnow;
axis off;
exportgraphics(gcf,'Figures/particle_n_128.jpeg','Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 2^10 particles
%% Number of particles and particles grid (their ratio must be a multiplee of 2)
r = 2^3; n = q*r; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Kerenel's FFT in both mean field and particle system
wHatXi = fft(wFun(xi));

%% Function handle for the particle system
sigma = sqrt(2*D);
a = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);
b = @(t,u) sigma*eye(length(u));

%% Function handle for the mean field
rhs = @(t,u) NeuralFieldFFT(u,wHatX,S,dx);

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
tspan = [0 35];
dt = 0.1; 
u0 = u0Fun(xi);
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

[XI,T] = meshgrid(xi,t);
surf(XI,T,U); shading interp; caxis(cax); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; % colorbar;
set(gcf,'position',pos);
drawnow;
axis off;
exportgraphics(gcf,'Figures/particle_n_1024.jpeg','Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 2^12 particles
%% Number of particles and particles grid (their ratio must be a multiplee of 2)
r = 2^6; n = q*r; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Kerenel's FFT in both mean field and particle system
wHatXi = fft(wFun(xi));

%% Function handle for the particle system
sigma = sqrt(2*D);
a = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);
b = @(t,u) sigma*eye(length(u));

%% Function handle for the mean field
rhs = @(t,u) NeuralFieldFFT(u,wHatX,S,dx);

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
tspan = [0 35];
dt = 0.1; 
u0 = u0Fun(xi);
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

[XI,T] = meshgrid(xi,t);
surf(XI,T,U); shading interp; caxis(cax); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; % colorbar;
set(gcf,'position',pos);
drawnow;
axis off;
exportgraphics(gcf,'Figures/particle_n_8129.jpeg','Resolution',300);



% rVals = 2.^[0:16];
% kVals = [0:12];

% %% Error plot
% eVals = zeros(length(rVals),length(kVals));
% uX = Umf(end,:)'; 
% for i = 1:length(rVals)

%   %% Number of particles and particles grid (their ratio must be a multiple of 2)
%   r = rVals(i); n = q*r; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%   %% Kerenel's FFT in both mean field and particle system
%   wHatXi = fft(wFun(xi));

%   %% Function handle for the particle system
%   sigma = sqrt(2*D);
%   a = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);
%   b = @(t,u) sigma*speye(length(u));

%   %% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
%   tspan = [0 35];
%   dt = 0.1; 
%   u0 = u0Fun(xi);
%   disp('TimeStep...');
%   [t,U] = EulerMaruyama(a,b,tspan,u0,dt);
%   disp('Done');

%   % subplot(3,2,[1,3]);
%   % [X,T] = meshgrid(x,tmf);
%   % surf(X,T,Umf); shading interp; caxis([-2 2]); view([0 90]); 
%   % shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
%   % axis tight; colorbar;
%   % drawnow;

%   % subplot(3,2,[2,4]);
%   % [XI,T] = meshgrid(xi,t);
%   % surf(XI,T,U); shading interp; caxis([-2 2]); view([0 90]); 
%   % shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
%   % axis tight; colorbar;

%   % %% Fourier coefficients
%   % uXi = U(end,:)'; 
%   uXi = U'; 
%   for l = 1:length(kVals)
%     k = kVals(l);
%     eVals(i,l) = abs( dxi*sum(exp(1i*k*pi/L*xi).*uXi) - dx*sum(exp(1i*k*pi/L*x).*uX) )
%   end

%   % subplot(3,2,[5,6]);
%   plot(x,uX,xi,uXi);
%   drawnow;
%   % pause

% end
% figure;
% nVals = q*rVals;
% loglog(nVals,eVals,'-*',nVals,nVals.^(-0.5));
% % legend({'k = 0',...
% %         'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5', 'k = 6', ... 
% %         'k = 7', 'k = 8', 'k = 9', 'k =10', 'k =11', 'k =12'});



