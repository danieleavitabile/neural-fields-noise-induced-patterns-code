
% % % % F(u) = -u + M f(u)

% M_ij = w(|xi -xj|)*2L/n
% K_ij = 
% M_ij = K_ij*2L/n

% % % p+(xi,xj), p-(xi,xj), p0(xi,xj) = 1-...
% % % K_ij in {-1, 0, 1}

% prob(K_ij = 1) =  p+(xi,xj)
% prob(K_ij =-1) =  p-(xi,xj)
% w(|xi -xj|) = p+(xi-xj)- p-(xi-xj)

%% Clean
clear all, close all, clc;

%% Reproducible random simulation
rng(1,"twister");

%% Parameters 
alpha = 10; theta  = 0.9; D = 0.5; B = 0.4; A = 1.0; u0Max = 5; u0Alpha = 0.25; L = 10*pi; % Turing

%% Funcion handle for the synaptic kernel
wFun = @(x) A*exp(-B*abs(x)).*(B*sin(abs(x)) + cos(x) );

%% Number of particles and particles grid (their ratio must be a multiplee of 2)
n = 2^12; dxi = 2*L/n; xi = -L+[0:n-1]'*dxi;

%% Deterministic matrix 
wVec = wFun(xi);

%% Define the positive and negative part of wFun, used later for probabilities
idPlus = find(wVec > 0);
idMinus = find(wVec < 0);
wPlus = 0*wVec; wPlus(idPlus) = wVec(idPlus);
wMinus = 0*wVec; wMinus(idMinus) = wVec(idMinus);

figure;
subplot(1,3,1);plot(xi, wVec);
subplot(1,3,2);plot(xi, wPlus);
subplot(1,3,3);plot(xi, wMinus);

%% Random extraction. Here I am assuming that the probability of an entry being 0 is 0.
%% Things will have to change for the case when p0 is not equal to 0.
wRand = zeros(size(xi));
for i = 1:n
  disp('*************')
  if ismember(i,idPlus)
    pPlus = wPlus(i);
    pMinus = 0;
    pZero = 1-pPlus-pMinus;
  else
    pPlus = 0;
    pMinus = -wMinus(i);
    pZero = 1-pPlus-pMinus;
  end
  % pPlus = wPlus(i)/(wPlus(i) - wMinus(i))
  % pMinus = -wMinus(i)/(wPlus(i) - wMinus(i))
  % pZero = 1 - pPlus -pMinus
  wRand(i) = randsample( [-1 0 1], 1, true, [pMinus pZero pPlus] );
end

figure, hold on;
plot(xi, wVec);
plot(xi, wRand * dxi, '*');
% plot(xi, wPlus, '*');
% plot(xi, wMinus, '*');

