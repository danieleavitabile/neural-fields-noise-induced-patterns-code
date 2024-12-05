function uAvg = ComputeAverage(u,xi,x)
   
  uAvg = zeros(size(x));
  q = length(x);
  dx = x(2)-x(1);
  n = length(xi);
  r = n/q;
  for i = 2:q-1
    uAvg(i)= sum(u((i-1)*r:(i+1)*r))/(2*dx);
  end

end
