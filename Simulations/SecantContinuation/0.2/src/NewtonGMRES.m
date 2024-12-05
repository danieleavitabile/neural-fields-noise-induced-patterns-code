function [x,resHist,flag] = NewtonGMRES(fhandle,jhandle,x0,options)

  % Rename parameters
  nltol       = options.nonlinTol;
  nlmaxit     = options.nonlinMaxIter;
  ltol        = options.linTol;
  restart     = options.linrestart;
  lmaxit      = options.linmaxit;
  display     = options.display;

  % Initialise iterations
  x = x0;
  f = feval(fhandle,x);
  neval = 1;
  res = norm(f,2);
  resHist = res;
  it = 0;

  % Displaying results
  if display
    DisplayIteration(it,neval,res);
  end

  % Main loop
  while ( (res > nltol) && (it < nlmaxit) )

    % Linear solve (gmres)
    [d,flag,relres,iter,resvec] = gmres(@(v) jhandle(x,v),-f,restart,ltol,lmaxit);
    if flag > 0
      warning(['GMRES did not converge. GMRES output flag = ' num2str(flag)]);
    end

    % Update solution
    x = x + d;

    % Book-keeping
    f = feval(fhandle, x); neval = neval + 1;
    res = norm(f,2); resHit = [resHist; res];
    it = it + 1;

    % Displaying results
    if display
      DisplayIteration(it,neval,res);
    end

  end

  % Output flag
  if res > nltol || isnan(res)
    flag = 0;
  else
    flag = 1;
  end

  function DisplayIteration(i,funceval,residual)
    if i == 0
      message = sprintf(['\n Start Newton-GMRES Iterations \n' ...
                           '   Iterations      Func-count      f(x)\n']);
      fprintf(message);
    end
    message = sprintf('%9d %16d %14.4e\n', i, funceval, residual);
    fprintf(message);
  end

end
  
