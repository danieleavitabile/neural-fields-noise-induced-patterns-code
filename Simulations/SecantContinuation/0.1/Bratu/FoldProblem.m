function [F,J] = FoldProblem(u,p)

  F = [p(1) - u(1)^2; p(2)-u(2)];

  if nargout > 1
    J = [-2*u(1) 0; 0 -1];
  end

end
