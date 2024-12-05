function W = GenerateSynapticMatrix(nodes,wFun,thresh,dataPath)

  %% Determine if a files must be saved
  if nargin < 4 || isempty(dataPath)
    saveFiles = false;
  else
    saveFiles = true;
  end

  % Get coordinates
  x = nodes(:,1); y = nodes(:,2); z = nodes(:,3); n = size(nodes,1);

  % % Generating synaptic matrix (sparsified with threshold)
  % % We are using here Euclidean distance, which is not correct on folded 
  % % surfaces (hence should change for real cortices)
  % D = sqrt((x - x').^2 + (y - y').^2 + (z - z').^2);
  % W = wFun(D);
  
  % % If a threshold is provided, then only elements above threshold
  % % in absolute value are effective, and a sparse matrix is given
  % innz = find(abs(W) > thresh);
  % W(setdiff(1:n^2,innz)) = 0;
  % sparse(W);

  % %% The code below could be vectorised (see above) but it is not viable when the
  % % problem is large, as D can not be stored in memory. The code below immediately 
  % % enforces sparsity, but it relies on for loops
  disp('Generate Synatpic Matrix...');
  tempnnz = 10000;
  rows = zeros(tempnnz);
  cols = zeros(size(rows));
  vals = zeros(size(rows));
  nnz = 0;
  for i = 1:n
    if mod(i,100) == 0
      disp(['Percentage of rows filled ' num2str(100*i/n)]);
    end
    for j = 1:n
      d = sqrt((x(i) - x(j)).^2 + (y(i) - y(j)).^2 + (z(i) - z(j)).^2);
      rho = wFun(d);
      if abs(rho) > thresh
        nnz = nnz + 1;
        rows(nnz) = i;
        cols(nnz) = j;
        vals(nnz) = rho;
      end
    end
  end
  idnnz = find(rows); 
  rows = rows(idnnz); 
  cols = cols(idnnz); 
  vals = vals(idnnz);
  W = sparse(rows,cols,vals,n,n);

  %% Save
  if saveFiles
    
    % Set filenames
    matFile = fullfile(dataPath,'synaptic-matrix.mat');
    datFile = fullfile(dataPath,'synaptic-matrix.dat');

    %% Save .mat file
    save(matFile,'W');
    disp(sprintf('Saved file %s',matFile));

    %% Save .dat file
    [row,col,wVal] = find(W);
    fileID = fopen(datFile,'w');
    fprintf(fileID,'#%12s %12s\n', 'numRows', 'numCols');
    fprintf(fileID,' %12i %12i\n', n, n);
    fprintf(fileID,'# Nonzero elements\n');
    fprintf(fileID,'#%12s %12s %20.12s\n', 'i', 'j', 'value');
    fprintf(fileID,' %12i %12i %20.12e\n', row, col, wVal);
    fclose(fileID);
    disp(sprintf('Saved file %s',datFile));

  end

end
