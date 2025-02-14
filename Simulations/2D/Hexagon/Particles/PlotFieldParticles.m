function fh = PlotFieldParticles(u,nodes,elem,opts)

  % Default options
  if ~isfield(opts,'xlim')
    opts.xlim = 'auto';
  end
  if ~isfield(opts,'ylim')
    opts.ylim = 'auto';
  end
  if ~isfield(opts,'zlim')
    opts.zlim = 'auto';
  end
  if ~isfield(opts,'clim')
    opts.clim = 'auto';
  end
  if ~isfield(opts,'view')
    opts.view = 3;
  end
  if ~isfield(opts,'axis')
    opts.axis = 'equal';
  end
  if ~isfield(opts,'markerSize')
    opts.markerSize = 50;
  end

  % Plot
  % fh = trisurf(elem,nodes(:,1),nodes(:,2),nodes(:,3),u);
  % fh.EdgeColor = '#4C566A'; %'none';
  % fh.EdgeAlpha = 0.4; %'none';
  e = ones(size(nodes(:,1)));
  fh = scatter(nodes(:,1), nodes(:,2),opts.markerSize*e,u,'filled');
  colorbar; 
  clim(opts.clim);
  view(opts.view); 
  axis(opts.axis)

  xlabel('x');
  ylabel('y');
  zlabel('z');

  xlim(opts.xlim);
  ylim(opts.ylim);
  zlim(opts.zlim);

end
