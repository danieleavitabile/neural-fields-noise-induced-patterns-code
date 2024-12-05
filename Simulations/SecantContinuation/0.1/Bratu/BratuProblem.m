function [F,J]=BratuProblem(u,p,L)

%% Rename parameters
alpha  = p(1);
lambda = p(2);
n2 = size(u,1);

%% Right-hand side
F = -alpha*L*u - lambda*exp(u);

%% Jacobian action
if nargout > 1
  J = -alpha*L - spdiags(lambda*exp(u),0,n2,n2);
end

end
