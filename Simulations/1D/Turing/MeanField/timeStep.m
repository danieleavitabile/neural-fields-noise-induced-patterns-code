%% Clean
clear all, close all, clc;

%% Parameters
% alpha = 10; theta  = 0.4; D = 10; B = 1.5; A = 30;
% alpha = 10; theta  = 0.4; D = 20; B = 1.5; A = 30;

% Increase noise to get patterns
% alpha = 10; theta  = 0.4; D = 0.02; B = 1.5; A = 7; % no patterns
alpha = 10; theta  = 0.4; D = 0.07; B = 1.5; A = 7;   % patterns

%% Funcion handles for the synaptic kernel
wFun = @(x) 1/sqrt(pi)*exp(-x.^2)- 1/(sqrt(pi)*B)*exp(-(x/B).^2); 
wHatFun = @(xi) ( exp(-xi.^2/4) - exp(-(xi*B).^2/4) );

%% Funcion handles for the firing rate function
phi = @(x) 0.5*(1+erf(x/sqrt(2)));
rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
S  = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
dS = @(u) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

%% Plot firing rate and kernel
figure;
x = linspace(-10,10,1000); u = linspace(0,5,1000);
subplot(3,3,1); 
plot(x,wFun(x)); xlabel('x'); legend({'w(x)'})
subplot(3,3,4); 
plot(u,S(u),u,dS(u)); xlabel('u'); legend({'S(u)','S''(u)'});

%% Plot Fourier Transform
xic = sqrt(8*log(B)/(B^2-1));
wHatc = wHatFun(xic);
Ac = 1/(wHatc*dS(0));

% Plot kernel fourier Transform. 
subplot(3,3,7);
xi = linspace(-6,6,1000);
plot(xi,wHatFun(xi),...
     xi, 1/(A*dS(0))*ones(size(xi)),...
     [-xic xic],wHatFun([-xic xic]),'*'...
); 

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
% u0 = 0.3*rand(size(x));
u0 = 0.3*cos(xic*x);
tspan = [0 100];
[t,U] = ode45(rhs,tspan,u0);

%% Plot state
subplot(3,3,[2 3 5 6 8 9]);
[X,T]= meshgrid(x,t);  
surf(X,T,U); view([0 90]); shading interp; caxis([-2 2]);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; colorbar;

%% Plot state on its own
figure;
[X,T]= meshgrid(x,t);  
surf(X,T,U); view([0 90]); shading interp; caxis([-2 2]);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
axis tight; colorbar;

p = D;
u = U(end,:)';
save('finalSate.mat','p','u');

