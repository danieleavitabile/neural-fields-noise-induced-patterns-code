%% Clean
clear all, close all, clc;

%% Parameters 
alpha = 10; theta  = 0.9; D = 0.1; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.25; L = 10*pi;

%% Funcion handle for the synaptic kernel
wFun = @(x) A*exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );

%% Funcion handles for the firing rates function (particles and mean field)
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
f = @(u) phi(alpha*(u-theta));
S = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));

%% Find homogeneous steady state
W0 = integral(wFun, -inf,inf);
G = @(u) -u +W0*S(u);
u = linspace(-1,1,1000);
plot(u,G(u));
uStar = 0.2;
uStar = fsolve(G,uStar)

%% Function handle for initial condition
u0Fun = @(x) uStar + u0Max./cosh(u0Alpha*x);

%% Mean field grid
q = 2^7; dx = 2*L/q; x = -L+[0:q-1]'*dx;

%% Number of particles and particles grid (their ratio must be a multiple of 2)
r = 2^3; n = q*r; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

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
figure;
subplot(2,2,[1,3]);
[X,T] = meshgrid(x,tmf);
surf(X,T,Umf); shading interp; caxis([-2 2]); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; colorbar;
drawnow;

%% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
tspan = [0 35];
dt = 0.1; 
u0 = u0Fun(xi);
[~,~,t,U] = EulerMaruyama(a,b,tspan,u0,dt);

subplot(2,2,[2,4]);
[XI,T] = meshgrid(xi,t);
surf(XI,T,U); shading interp; caxis([-2 2]); view([0 90]); 
shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; colorbar;

rVals = 2.^[0:14];
kVals = [0:20];

%% Error plot
eVals = zeros(length(rVals),length(kVals));
uX = Umf(end,:)'; 
for i = 1:length(rVals)

  %% Number of particles and particles grid (their ratio must be a multiple of 2)
  r = rVals(i); n = q*r; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

  %% Kerenel's FFT in both mean field and particle system
  wHatXi = fft(wFun(xi));

  %% Function handle for the particle system
  sigma = sqrt(2*D);
  a = @(t,u) NeuralFieldFFT(u,wHatXi,f,dxi);
  b = @(t,u) sigma*speye(length(u));

  %% Euler-Maruyama function handles for model dU = a(t,U)dt + b(t,U)dW
  tspan = [0 35];
  dt = 0.1; 
  u0 = u0Fun(xi);
  disp('TimeStep...');
  [t,U] = EulerMaruyama(a,b,tspan,u0,dt);
  disp('Done');

  % subplot(3,2,[1,3]);
  % [X,T] = meshgrid(x,tmf);
  % surf(X,T,Umf); shading interp; caxis([-2 2]); view([0 90]); 
  % shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
  % axis tight; colorbar;
  % drawnow;

  % subplot(3,2,[2,4]);
  % [XI,T] = meshgrid(xi,t);
  % surf(XI,T,U); shading interp; caxis([-2 2]); view([0 90]); 
  % shading interp; xlabel('x'); ylabel('t'); zlabel('u(x,t)');
  % axis tight; colorbar;

  % %% Fourier coefficients
  % uXi = U(end,:)'; 
  uXi = U'; 
  for l = 1:length(kVals)
    k = kVals(l);
    eVals(i,l) = abs( dxi*sum(exp(1i*k*pi/L*xi).*uXi) - dx*sum(exp(1i*k*pi/L*x).*uX) )
  end

  % subplot(3,2,[5,6]);
  figure, hold on;
  % title(['n =' num2str(n)]);
  greyNord = '#D8DEE9';
  blueNord = '#5E81AC';
  redNord = '#B48EAD';
  greenNord = '#A3BE8C';

  plot(xi,uXi,'color',greenNord,'linewidth',1);
  plot(x,uX,'color',blueNord,'linewidth',3);
  box on; axis tight; xlim = ([-L L]); ylim = ([-2.5 4]);
  fileName = sprintf('Figures/profile_n_%07i.eps',n);
  set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
  % exportgraphics(gcf,fileName,'resolution',600);
  set(gcf, 'color', 'none');    
  set(gca, 'color', 'none');
  exportgraphics(gcf,fileName,'ContentType','vector','BackgroundColor','none');
  drawnow;
  % pause

end
figure;
nVals = q*rVals; axis tight; ylim = ([1e-4 1e2]);
loglog(nVals,eVals,'*-',nVals,nVals.^(-0.5));
% legend({'k = 0',...
%         'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5', 'k = 6', ... 
%         'k = 7', 'k = 8', 'k = 9', 'k =10', 'k =11', 'k =12'});

