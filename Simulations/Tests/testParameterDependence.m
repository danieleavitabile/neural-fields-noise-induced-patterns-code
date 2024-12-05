%% Clean
clear all, close all, clc;

%% Parameters
% alpha = 10; theta  = 0.4; D = 10; B = 1.5; A = 30;
alpha = 10; theta  = 0.4; D = 20; B = 1.5; A = 30;

%% Funcion handles for the synaptic kernel
wFun = @(x) 1/sqrt(pi)*exp(-x.^2)- 1/(sqrt(pi)*B)*exp(-(x/B).^2); 
wHatFun = @(xi) ( exp(-xi.^2/4) - exp(-(xi*B).^2/4) );

%% Funcion handles for the firing rate function
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);

%% Dependence on alpha
S  = @(u,D,alpha,theta) phi(alpha*(u-theta)./sqrt(1+alpha.^2*D));
dS = @(u,D,alpha,theta) rho(alpha*(u-theta)./sqrt(1+alpha.^2*D),0,1)*alpha./sqrt(1+alpha.^2*D);
D = linspace(0,20,10000);
figure, hold on;
plot(D,dS(0,D,0.1*alpha,theta));
plot(D,dS(0,D,alpha,theta));
plot(D,dS(0,D,10*alpha,theta));
plot(D,dS(0,D,100*alpha,theta));
xlabel('D');
ylabel('F''(0)');
legend({'\alpha = 1', '\alpha = 10', '\alpha = 100', '\alpha=1000'});

%% Dependence on theta
S  = @(u,D,alpha,theta) phi(alpha*(u-theta)./sqrt(1+alpha.^2*D));
dS = @(u,D,alpha,theta) rho(alpha*(u-theta)./sqrt(1+alpha.^2*D),0,1)*alpha./sqrt(1+alpha.^2*D);
theta = linspace(0,10,10000);
figure, hold on;
plot(D,dS(0,0.1,alpha,theta));
plot(D,dS(0,2.5,alpha,theta));
plot(D,dS(0,5,alpha,theta));
xlabel('\theta');
ylabel('F''(0)');
legend({'D = 0.1', 'D = 2.5', 'D = 5'});

