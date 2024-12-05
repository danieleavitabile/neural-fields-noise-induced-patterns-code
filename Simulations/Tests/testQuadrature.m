clear all, close all, clc;

u = linspace(-10,10,1000);
D = 10;
mu = 10;
h = 1;

f = @(x) 1./(1+exp(-mu*(x-h)));
rho = @(x,u,D) 1/(2*pi*D)^0.5 * exp(-(x-u).^2/(2*D));

I1 = zeros(size(u));
I2 = zeros(size(u));
for k = 1:length(u)
  g1 = @(x) f(x).*rho(x,u(k),D);
  g2 = @(x) f(x).*rho(x,u(k),10*D);
  I1(k) = integral(g1,-Inf,Inf);
  I2(k) = integral(g2,-Inf,Inf);
end

% figure; plot(u, f(u));
% figure; plot(u, rho(u,1));
% figure; plot(u, f(u), u,I1,u,I2);

[nodesGH,weightsGH] = GaussHermite(120);

I3 = zeros(size(u));
for k = 1:length(u)
  I3(k)= MeanRate(f,u(k),D,nodesGH,weightsGH);
end

plot(u, I1, u, I3,'*')

% u = 10;
% integral(@(x) f(x).*rho(x,u,D),-Inf,Inf)
% MeanRate(f,u,D,nodesGH, weightsGH)

% MeanRate(@(x) cos(x), 0,1/2,weightsGH, nodesGH)



% clear all, close all, clc;
% n   = 7;
% [x,w] = gauss_her(n);
% I   = dot(w,cos(x))
% [x,w] = GaussHermite(7);
% I   = MeanRate(@cos,0,1/2,x,w)


% function [x,w] = gauss_her(n);
% % Build Hn (Hermite Poly of order n) by recurrence :
% h(1,1)=1;
% h(2,1:2)=[2 0];
% for k=2:n
%     h(k+1,1:k+1)=2*[h(k,1:k) 0]-2*(k-1)*[0 0 h(k-1,1:k-1)];
% end


% Hn       = h(n+1,:);
% Hn_deriv = polyder(Hn);

% x        = roots(Hn);
% w        = 2^(n+1)*factorial(n)*sqrt(pi)  ./...
%            (polyval(Hn_deriv,x)).^2;
% end
