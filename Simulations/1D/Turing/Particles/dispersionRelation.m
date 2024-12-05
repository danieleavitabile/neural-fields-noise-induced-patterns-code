%% Clean
clear all, close all, clc;

%% Parameters 
% alpha = 10; theta  = 0.4; D = 0.01; B = 1.5; A = 7;   
alpha = 10; theta  = 0.4; D = 0.065; B = 1.5; A = 7;   % patterns

%% Funcion handles for the synaptic kernel
wFun = @(x) A*(1/sqrt(pi)*exp(-x.^2)- 1/(sqrt(pi)*B)*exp(-(x/B).^2)); 

%% Firing rate
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
f = @(u) phi(alpha*(u-theta));
S  = @(u,D) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u,D) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Number of particles and spatial grid
n = 2^12; L = 10*pi; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Homogeneous steady state function
a0 = A*( erf(L)-erf(L/B) );
m0Fun = @(m) -m + S(m,D)*a0;
m0 = fsolve(m0Fun,0);

figure, hold on;
kVals = 0:50;
% AKVals = zeros(size(kVals));
lambdaKVals = zeros(size(kVals));
for j = 1:length(kVals)
  Ak = integral(@(x) (2*L)^(-1)*wFun(x).*exp(1i*kVals(j)*pi*x/L),-L,L);
  lambdaKVals(j) = -1 + dS(m0,D)*2*L*Ak;
end
plot(kVals,real(lambdaKVals),'.','DisplayName','D = 0.065');
plot(kVals,0*kVals);
xlabel('k'); ylabel('Re lk');
drawnow;

alpha = 10; theta  = 0.4; D = 0.060; B = 1.5; A = 7;   % patterns

%% Funcion handles for the synaptic kernel
wFun = @(x) A*(1/sqrt(pi)*exp(-x.^2)- 1/(sqrt(pi)*B)*exp(-(x/B).^2)); 

%% Firing rate
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
f = @(u) phi(alpha*(u-theta));
S  = @(u,D) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u,D) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Homogeneous steady state function
a0 = A*( erf(L)-erf(L/B) );
m0Fun = @(m) -m + S(m,D)*a0;
m0 = fsolve(m0Fun,0);

for j = 1:length(kVals)
  Ak = integral(@(x) (2*L)^(-1)*wFun(x).*exp(1i*kVals(j)*pi*x/L),-L,L);
  lambdaKVals(j) = -1 + dS(m0,D)*2*L*Ak;
end
plot(kVals,real(lambdaKVals),'.','DisplayName','D = 0.060');
xlabel('k'); ylabel('Re lk');
legend;
drawnow;

DVals = 0:0.01:0.5;
kVals = 0:20;
lambdaVals = zeros(length(DVals),length(kVals));
for i = 1:length(DVals)

  D = DVals(i);
  m0Fun = @(m) -m + S(m,D)*a0;
  m0 = fsolve(m0Fun,0);
  for j = 1:length(kVals)
    Ak = integral(@(x) (2*L)^(-1)*wFun(x).*exp(1i*kVals(j)*pi*x/L),-L,L);
    lambdaKVals(i,j) = real(-1 + dS(m0,D)*2*L*Ak);
  end

end

figure, hold on;
plot(DVals,0*DVals);
for j = 1:length(kVals)
  grey = [.5 .5 .5];
  lw = 0.01;
  plot(DVals,lambdaKVals(:,j),'color', grey, 'linewidth', lw);
end
xlabel('sigma');
ylabel('Re l_k');
