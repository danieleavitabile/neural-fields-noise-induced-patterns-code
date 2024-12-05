function fh = PlotStateParticles(z,nodes,elem,opts)

  n = size(z,1)/2; iU = 1:n; iV = n+iU;
  u = z(iU); v = z(iV);

  fh = gcf;
  subplot(1,2,1);
  PlotFieldParticles(u,nodes,elem,opts.U);

  subplot(1,2,2);
  PlotFieldParticles(v,nodes,elem,opts.V);

end
