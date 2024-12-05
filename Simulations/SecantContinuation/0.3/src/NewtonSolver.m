function [x,res,flag,it] = NewtonSolver(fhandle,x0,options)

  % Rename parameters
  nltol       = options.nonlinTol;
  nlmaxit     = options.nonlinMaxIter;
  display     = options.display;

  % Initialise iterations
  x = x0;
  [f,J] = feval(fhandle,x);
  neval = 1;
  res = norm(f,2);
  itHist = [res neval];
  it = 0;

  % Displaying results
  if display
    DisplayIteration(it,neval,res);
  end

  % Main loop
  while ( (res > nltol) && (it < nlmaxit) )

    % Linear solve (direct)
    %J = Jacobian(x);
    tic
    disp('Linear Solve...');
    size(x), size(f), size(J)
    d = J\(-f);
    toc

    % Update solution
    x = x + d;

    % Book-keeping
    [f,J] = feval(fhandle, x); neval = neval + 1;
    res = norm(f,2);
    itHist = [itHist; [res neval]];
    it = it + 1;

    % Displaying results
    if display
      DisplayIteration(it,neval,res);
    end

  end

  % Output flag
  if res > nltol || isnan(res)
    flag = -1;
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

