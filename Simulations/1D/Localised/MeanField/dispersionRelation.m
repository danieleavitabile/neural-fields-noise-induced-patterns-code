%% Clean
clear all, close all, clc;

%% Load homogeneous steady state (middle point, which is used as initial guess for Newton)
sol = load('Homogeneous/solution_0000020.mat');
m0 = sol.u(length(sol.u)/2);

%% Parameters lifted from the converged solution (and others completed/overridden)
alpha = sol.p(1); theta = sol.p(2); D = sol.p(3); A = sol.p(4); B = 0.4; u0Max = 5; u0Alpha = 0.25; L = 10*pi;

%% Override D
D = 0.45; sigma = sqrt(2*D)

%% Funcion handles for the synaptic kernel
wFun = @(x) A*exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );

%% Firing rate
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
f = @(u) phi(alpha*(u-theta));
S  = @(u,D) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u,D) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Number of particles and spatial grid
n = 2^12; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Converge initial guess coming from spatially-extended system to steady state 
a0 = integral(@(x) wFun(x),-L,L);
m0Fun = @(m) -m + S(m,D)*a0;
m0
m0 = fsolve(m0Fun,m0)

figure, hold on;
kVals = 0:50;
% AKVals = zeros(size(kVals));
lambdaKVals = zeros(size(kVals));
for j = 1:length(kVals)
  Ak = integral(@(x) (2*L)^(-1)*wFun(x).*exp(1i*kVals(j)*pi*x/L),-L,L);
  lambdaKVals(j) = -1 + dS(m0,D)*2*L*Ak;
end
plot(kVals,real(lambdaKVals),'.','DisplayName',['D = ', num2str(D)]);
plot(kVals,0*kVals);
xlabel('k'); ylabel('Re lk');
drawnow;

%% Override D
D = 0.41; sigma = sqrt(2*D)

%% Funcion handles for the synaptic kernel
wFun = @(x) A*exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );

%% Firing rate
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
f = @(u) phi(alpha*(u-theta));
S  = @(u,D) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u,D) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Number of particles and spatial grid
n = 2^12; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Converge initial guess coming from spatially-extended system to steady state 
a0 = integral(@(x) wFun(x),-L,L);
m0Fun = @(m) -m + S(m,D)*a0;
m0
m0 = fsolve(m0Fun,m0)

kVals = 0:50;
lambdaKVals = zeros(size(kVals));
for j = 1:length(kVals)
  Ak = integral(@(x) (2*L)^(-1)*wFun(x).*exp(1i*kVals(j)*pi*x/L),-L,L);
  lambdaKVals(j) = -1 + dS(m0,D)*2*L*Ak;
end
plot(kVals,real(lambdaKVals),'.','DisplayName',['D = ', num2str(D)]);
xlabel('k'); ylabel('Re lk');
legend;
drawnow;

% % alpha = 10; theta  = 0.4; D = 0.060; B = 1.5; A = 7;   % patterns

% % %% Funcion handles for the synaptic kernel
% % wFun = @(x) A*(1/sqrt(pi)*exp(-x.^2)- 1/(sqrt(pi)*B)*exp(-(x/B).^2)); 

% % %% Firing rate
% % phi = @(x) 0.5*(1+erf(x/sqrt(2)));
% % rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
% % f = @(u) phi(alpha*(u-theta));
% % S  = @(u,D) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
% % dS = @(u,D) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

% % %% Homogeneous steady state function
% % a0 = A*( erf(L)-erf(L/B) );
% % m0Fun = @(m) -m + S(m,D)*a0;
% % m0 = fsolve(m0Fun,0);

% % for j = 1:length(kVals)
% %   Ak = integral(@(x) (2*L)^(-1)*wFun(x).*exp(1i*kVals(j)*pi*x/L),-L,L);
% %   lambdaKVals(j) = -1 + dS(m0,D)*2*L*Ak;
% % end
% % plot(kVals,real(lambdaKVals),'.','DisplayName','D = 0.060');
% % xlabel('k'); ylabel('Re lk');
% % legend;
% % drawnow;

% % DVals = 0:0.01:0.5;
% % kVals = 0:20;
% % lambdaVals = zeros(length(DVals),length(kVals));
% % for i = 1:length(DVals)

% %   D = DVals(i);
% %   m0Fun = @(m) -m + S(m,D)*a0;
% %   m0 = fsolve(m0Fun,0);
% %   for j = 1:length(kVals)
% %     Ak = integral(@(x) (2*L)^(-1)*wFun(x).*exp(1i*kVals(j)*pi*x/L),-L,L);
% %     lambdaKVals(i,j) = real(-1 + dS(m0,D)*2*L*Ak);
% %   end

% % end

% % figure, hold on;
% % plot(DVals,0*DVals);
% % for j = 1:length(kVals)
% %   grey = [.5 .5 .5];
% %   lw = 0.01;
% %   plot(DVals,lambdaKVals(:,j),'color', grey, 'linewidth', lw);
% % end
% % xlabel('sigma');
% % ylabel('Re l_k');
