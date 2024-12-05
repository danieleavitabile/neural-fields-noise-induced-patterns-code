% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = './';
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'CMU Bright';
opts.fontSize   = 10;

% create new figure
fig = gcf

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;

% set text properties
set(fig.Children, ...
    'FontName',    opts.fontType, ...
    'FontSize',    opts.fontSize);

% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

% export to png
fig.PaperPositionMode   = 'auto';
fileName = [opts.saveFolder 'my_figure.eps'];
% print([opts.saveFolder 'my_figure'], '-dpdf')
exportgraphics(fig, fileName, 'ContentType', 'vector');
