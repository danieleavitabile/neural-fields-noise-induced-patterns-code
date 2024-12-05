function plotHandle = PlotSolution(u,p,parentHandle,X,Y,view2D)

    if isempty(parentHandle)
      scrsz = get(0,'ScreenSize');
      plotHandle = figure('Position',[3/4*scrsz(3) scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
      parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
    end

    [ny,nx] = size(X);
    U = reshape(u,ny,nx);

    figure(parentHandle);
    surf(X,Y,U);
    titleString = 'U';
    title(titleString);
    shading interp;
    if view2D
      view([0 90]); axis tight;
    end
    caxis([0 8]);
    axis square;

    print -dtiff state.tiff

end

