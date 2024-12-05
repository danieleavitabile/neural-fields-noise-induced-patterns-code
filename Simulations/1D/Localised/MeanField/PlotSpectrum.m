function plotHandle = PlotSpectrum(d,p,parentHandle)

   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[3*scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end

   figure(parentHandle);
   cla, hold on;

   % reSpan = [-2 2];
   % imSpan = [-2 2];
   reSpan = [min(real(d)) max(real(d))];
   imSpan = [min(imag(d)) max(imag(d))];

   plot(real(d),imag(d),'.','MarkerSize',20);
   plot(zeros(1,100),linspace(imSpan(1),imSpan(2)));
   % xlim(reSpan); ylim(imSpan); 
   grid on;
   hold off;

   %% Save
   % print -dtiff spectrum.tiff

end
