clear all, close all, clc;

phi = @(x) 0.5*(1+erf(x/sqrt(2)));

rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);

u = linspace(-10,10,1000);
D = 10;
alpha = 10;
h = 1;

I = zeros(size(u));
for k = 1:length(u)
  g = @(x) phi(alpha*(x-h)).*rho(x,u(k),D);
  I(k) = integral(g,-Inf,Inf);
end

phiExp = @(u) phi(alpha*(u-h)./sqrt(1+alpha^2*D));

figure; plot(u, I,u, phiExp(u),'*');


