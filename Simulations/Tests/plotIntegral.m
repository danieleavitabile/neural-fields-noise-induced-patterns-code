clear all, close all, clc;

u = linspace(-10,10,1000);
D = 10;
mu = 10;
h = 1;

f = @(x) 1./(1+exp(-mu*(x-h)));
rho = @(x,u,D) 1/(2*pi*D)^0.5 * exp(-(x-u).^2/(2*D));

I = zeros(size(u));
for k = 1:length(u)
  g1 = @(x) f(x).*rho(x,u(k),D);
  g2 = @(x) f(x).*rho(x,u(k),10*D);
  I1(k) = integral(g1,-Inf,Inf);
  I2(k) = integral(g2,-Inf,Inf);
end

% figure; plot(u, f(u));
% figure; plot(u, rho(u,1));
figure; plot(u, f(u), u,I1,u,I2);

