%% Clean
clear all, close all, clc;

%% Parameters
% alpha = 10; theta  = 0.9; D = 0.1; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 4;
alpha = 10; theta  = 0.9; D = 0.1; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.18;
% alpha = 10; theta  = 0.9; D = 0.1; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.04;
% alpha = 10; theta  = 0.9; D = 0.1; B = 0.4; A = 1.0; u0Max = 3; u0Alpha = 0.4;

%% Funcion handles for the synaptic kernel
% wFun = @(x) 1/sqrt(pi)*exp(-x.^2)- 1/(sqrt(pi)*B)*exp(-(x/B).^2); 
% wHatFun = @(xi) ( exp(-xi.^2/4) - exp(-(xi*B).^2/4) );
wFun = @(x) exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );
wHatFun = @(xi) ...
           2*B*(B^2 -xi.^2 +1)./(B^4+2*B^2*xi.^2 + 2*B^2 + xi.^4 -2*xi.^2 +1)...
         + B./(B^2+(xi-1).^2) + B./(B^2+(xi+1).^2);

%% Funcion handles for the firing rate function
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
S  = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

W0 = integral(wFun, -inf,inf);
G = @(u) -u +W0*S(u);
u = linspace(-1,1,1000);
plot(u,G(u));
uStar = 0.2;
uStar = fsolve(G,uStar)

%% Plot firing rate and kernel
figure;
x = linspace(-30,30,1000); u = linspace(0,5,1000);
subplot(3,3,1); 
plot(x,wFun(x)); xlabel('x'); legend({'w(x)'})
subplot(3,3,4); 
plot(u,S(u),u,dS(u)); xlabel('u'); legend({'S(u)','S''(u)'});

%% Plot Fourier Transform
% xic = sqrt(8*log(B)/(B^2-1));
xi = linspace(-12,12,1000);
wHatc = max(wHatFun(xi));
Ac = 1/(wHatc*dS(uStar));

% % Plot kernel fourier Transform. 
subplot(3,3,7);
xi = linspace(-12,12,1000);
plot(xi,wHatFun(xi),...
     xi, 1/(A*dS(uStar))*ones(size(xi))); 

%% Spatial grid
n = 2^10; L = 10*pi; dx = 2*L/n; x = -L+[0:n-1]'*dx;

%% Form matrix (ring geometry)
M = zeros(n,n);
y = wFun(x)*dx;
iRows = 1:n;
iShift = -n/2:n/2-1;
for i = 1:n
  M(iRows(i),:) = circshift(y, iShift(i));
end

%% Form the neural field RHS
rhs = @(t,u) -u + A*M*S(u);

%% Initial condition, and time step
% u0 = uStar + 0.1*rand(size(x));
% u0 = 0.3*cos(xic*x);
u0 = uStar + u0Max./cosh(u0Alpha*x);
tspan = [0 100];
[t,U] = ode45(rhs,tspan,u0);

%% Plot state
subplot(3,3,[2 3 5 6 8 9]);
[X,T]= meshgrid(x,t);  
surf(X,T,U); view([0 90]); shading interp;
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
caxis([0 4]);
axis tight; colorbar;
