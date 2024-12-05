function [tFinal,uFinal,tHist,uHist,dwHist] = EulerMaruyama(a,b,tspan,u0,dt)

  % Preprocess input
  if nargout > 2
    saveHist = true;
  else
    saveHist = false;
  end
  if nargout > 4
    saveNoise = true;
  else
    saveNoise = false;
  end
  if isrow(u0) 
    error('u0 must be a column vector');
  end

  % Instantiate history vectors
  nt = (tspan(2)-tspan(1))/dt - 1;
  n = size(u0,1);
  if saveHist
    tHist = zeros(nt+1,1);
    uHist = zeros(nt+1,n);
  end
  if saveNoise then
    dwHist = zeros(nt,n);
  end

  % Initialise
  t = tspan(1); u = u0;
  if saveHist
    tHist(1) = t; uHist(1,:) = u';
  end

  % Main loop
  for i = 1:nt

    % Take one step
    dw = normrnd(0,sqrt(dt),n,1);
    u = u + a(t,u)*dt + b(t,u)*dw;
    t = t + dt;

    % Save history
    if saveHist
      tHist(i+1) = t;
      uHist(i+1,:) = u';
    end
    if saveNoise
      dwHist(i,:) = dw';
    end

  end

  tFinal = t;
  uFinal = u';

end
