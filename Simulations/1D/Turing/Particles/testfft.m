clear all; close all; clc;

L = 3*pi;

q = 2^7; dx = 2*L/q; x = -L+[0:q-1]'*dx;

yx = sin(pi/L*x) + 0.01*rand(size(x));

% n = 2^3*q; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;
% yxi = sin(pi/L*xi) + 0.1*rand(size(xi));

% figure
% subplot(2,2,[1 2]);
% plot(x,yx,xi,yxi);

% yxHat = fft(yx);

% yxHat = dx/q*yxHat(1:q/2+1); 
% yxPow = abs(yxHat).^2; yxPow(2:end-1) = 2*yxPow(2:end-1);
% omegax = 0:dx:L;

% subplot(2,2,[3 4]);
% plot(omegax, yxPow);

