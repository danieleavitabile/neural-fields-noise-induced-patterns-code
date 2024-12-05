function F = NeuralFieldFFT(u,wHat,rateFun,dx)

  %% Rename parameters
  n = size(u,1);

  %% Firing rate function
  fHat = fft(rateFun(u));

  %% Right-hand side
  F = zeros(n,1);
  F = - u + dx*ifftshift(real(ifft(fHat .* wHat)));
  

end
