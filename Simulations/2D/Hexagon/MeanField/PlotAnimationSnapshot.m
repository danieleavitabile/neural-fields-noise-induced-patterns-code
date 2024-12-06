function PlotAnimationSnapshot(id,t,uHist,nodes,elem,opts)

  % Time steps and number of nodes
  disp(['t =' num2str(t(id))] );

  % Plot initial state
  fh = PlotField(uHist(id,:)',nodes,elem,opts);
  set(gca,'visible','off');
  colorbar off;
  shading interp;

end
