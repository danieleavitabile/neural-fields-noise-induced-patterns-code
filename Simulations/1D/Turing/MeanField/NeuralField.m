function [F,JPsi] = NeuralField(u,p,wHat,dx,psi);

  % Rename parameters
  alpha  = p(1);
  theta  = p(2);
  D      = p(3);

  % Firing rate parameters
  phi = @(x) 0.5*(1+erf(x/sqrt(2)));
  rho = @(x,u,D) exp(-(x-u).^2./(2*D) )/sqrt(2*pi*D);
  S  = @(u) phi(alpha*(u-theta)./sqrt(1+alpha^2*D));
  dS = @(u) rho(alpha*(u-theta)./sqrt(1+alpha^2*D),0,1)*alpha/sqrt(1+alpha^2*D);

  % Convolution between w and r
  rConv = dx*ifftshift(real(ifft( fft(S(u)) .* wHat )));

  % Right-hand side
  F = -u +rConv;
  if any(~isreal(F))
    alpha
    theta
    D
    pause
  end

  % Jacobian action
  if nargout > 1 && ~isempty(psi)
    psiConv = dx*ifftshift(real(ifft( fft(dS(u).*psi) .* wHat )));
    JPsi = -psi + psiConv;
  if any(~isreal(JPsi))
    alpha
    theta
    D
    pause
  end
  end


end
