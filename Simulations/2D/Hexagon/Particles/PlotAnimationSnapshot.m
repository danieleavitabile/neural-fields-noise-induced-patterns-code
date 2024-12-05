function PlotAnimationSnapshot(id,t,uHist,nodes,elem,opts)

  % Time steps and number of nodes
  disp(['t =' num2str(t(id))] );

  % Plot initial state
  fh = PlotFieldParticles(uHist(id,:)',nodes,elem,opts);
  set(gca,'visible','off');
  colorbar off;
  shading interp;

end
