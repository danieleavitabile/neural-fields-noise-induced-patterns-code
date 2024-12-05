function F = SolutionMeasures(step,u,p,hx,P)

  %% Compute 2-norm
  l2Norm = sqrt(sum(hx * u.^2)/P);

  %% Allocate
  F = zeros(1,5);

  %% Assign branch variables
  F = [l2Norm max(u) min(u)];
  
end
