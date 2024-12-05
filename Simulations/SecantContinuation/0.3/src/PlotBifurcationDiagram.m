function plotHandle = PlotBifurcationDiagram(branchFile,idVar)

  if isempty(branchFile)
    [fileName,pathName] = uigetfile('*.mat','Select the branch file');
    branchFile = [pathName fileName];
  end

  plotHandle = figure; hold on;

  if ~iscell(branchFile)
    branchFile = {branchFile};
  end
  numBranches = length(branchFile);

  for b = 1:numBranches

    branchData = load(char(branchFile(b)));
    branch = branchData.branch;
    iUnstab = find(branch(:,2) > 0);
    iStab   = find(branch(:,2) == 0); 
    iNoStab = find(branch(:,2) == -1); 
    if ~isempty(iNoStab)
      plot(branch(:,3),branch(:,idVar),'k-');
    end
    if ~isempty(iUnstab)
      plot(branch(:,3),branch(:,idVar),'b-');
    end
    if ~isempty(iStab)
      jumps = find( diff(iStab) > 1 )
      if isempty(jumps)
	plot(branch(:,3),branch(:,idVar),'b-','LineWidth',1.5);
      else
	numJumps = length(jumps);
	numPatches = numJumps + 1;
	for p = 1:numPatches
	  if p == 1
	      pIni = iStab(1); pEnd = iStab(jumps(p));
	  elseif p == numPatches
	      pIni = iStab(jumps(p-1) + 1); pEnd = iStab(end);
	  else
		  pIni = iStab(jumps(p-1) + 1); pEnd = iStab(jumps(p));
	  end
	      plot(branch(pIni:pEnd,3),branch(pIni:pEnd,idVar),'b-','LineWidth',1.5);
	end
      end
    end

  end

  hold off;
  set(gca,'Box','on');

  print -depsc branch.eps;

end

