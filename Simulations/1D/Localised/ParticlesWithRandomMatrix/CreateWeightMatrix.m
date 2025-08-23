function WMat = CreateWeightMatrix(wFun,xi)

  %% Assumes wFun is a function in the interval [-1,1]
  %% Assumes n = length(xi) is even
  %% Assumes xi coves in [-L, L-dx], and a convolutional structure
  %% We use the convolutional structure of the matrix. For the numbering of rows and
  %% columns see page 3, Tutorial 1 of
  %% https://github.com/danieleavitabile/numerical-analysis-mathematical-neuroscience/
  %% 
  %% We define probabilities using the values w(0-xj) which corresponds to row n/2 of
  %% the matrix. The probabilities in all other rows are then just circularly shifted.
  %% We exploit this property by: iterating on the columns j, extract n samples
  %% (whhich will be used) in the whole column j, and circularly shifted. This saves
  %% considerable computational time, in Matlab, with respect to running a double for
  %% loop.

  % Rename parameters
  n = length(xi); dxi = xi(2) - xi(1); 

  % Store values of w
  wVec = wFun(xi);

  % Define the positive and negative part of w, used later for probabilities
  idPlus = find(wVec > 0);
  idMinus = find(wVec < 0);
  wPlus = 0*wVec; wPlus(idPlus) = wVec(idPlus);
  wMinus = 0*wVec; wMinus(idMinus) = wVec(idMinus);
  
  % For each column
  for j = 1:n

    % Determine probabilities
    wVal = wVec(j);
    if wVal > 0
      pPlus = wVal;
      pMinus = 0;
      pZero = 1-pPlus-pMinus;
    elseif wVal < 0
      pPlus = 0;
      pMinus = -wVal;
      pZero = 1-pPlus-pMinus;
    else
      pPlus = 0;
      pMinus = 0;
      pZero = 1;
    end

    % Extrant n samples for a column j (these will be appropriately shifted later)
    KMat(:,j) = randsample( [-1 0 1], n, true, [pMinus pZero pPlus] );

  end

  % The samples now correspond to w(0-xj) for all j. We scale them and circularly
  % shift, as appropriate for the convolutional structure of the problem
  for i = 1:n
    WMat(i,:) = circshift( KMat(i,:), -n/2-1+i )*dxi;
  end
  
  % Declare matrix as sparse
  WMat = sparse(WMat);

end
