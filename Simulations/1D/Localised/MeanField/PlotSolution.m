function plotHandle = PlotSolution(x,u,p,parentHandle);

   %% Position and eventually grab figure
   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end

   figure(parentHandle);

   plot(x,u,'.-');

   %% Save
   % print -depsc state.eps
   % print -dtiff state.tiff

end

