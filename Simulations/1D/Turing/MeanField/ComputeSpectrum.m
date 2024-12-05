function [V,D] = ComputeSpectrum(u,p,wHat,dx);

    %% Linear function handle
    jhandle = @(v) NeuralFieldJacobianAction(u,p,wHat,dx,v);

    %% Call eigenvalue solver
    eigsOpts = [];

    [V,D,flag] = eigs(jhandle,size(u,1),10,'lr');

    if flag ~= 0
      warning('Not all eigenvalues converged');
    end

end
