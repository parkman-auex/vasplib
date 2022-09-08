function varargout = WilsonLoopPlot(BFCAR,klist_l,options)
arguments
    BFCAR;
    klist_l;
    options.ycut = [-pi,pi];
    options.shift = 0;
    options.Color = [rand rand rand];
    options.LineSpec = 'o';
    options.title = 'WilsonLoop';
    options.xlabel='k_{envolution}';
    options.ylabel=["WannierCenter";"(BerryPhase1D)"];
    options.ax = handle([]);
end
if options.ycut(1) >= 0
    BFCAR(BFCAR(:,:)<0) = 2*pi+BFCAR(BFCAR(:,:)<0);
end
if isempty(options.ax)
    Fig = vasplib_plot.create_figure('Position',[0.2,0.2,0.6,0.6]);
    ax = Fig.axes(1);
else
    ax = options.ax;
end
if length(klist_l)==1
    BFCAR = [BFCAR,BFCAR];
    klist_l = [0,klist_l];
end
ax = vasplib_plot.bandplot(BFCAR+options.shift,options.ycut+options.shift,klist_l,[],[],...
    'LineSpec',options.LineSpec,...
    'Color',options.Color,...
    'xlabel',options.xlabel,...
    'ylabel',options.ylabel, ...
    'ax',ax,...
    'title',options.title  ...
    );
set(ax,'xtick',-2*pi:pi/2:2*pi,...
    'xticklabel',["-2\pi","-3\pi/2","-\pi","-\pi/2",...
    "0","\pi/2","\pi","3\pi/2","2\pi"]);
set(ax,'ytick',-2*pi:pi/2:2*pi,...
    'yticklabel',["-2\pi","-3\pi/2","-\pi","-\pi/2",...
    "0","\pi/2","\pi","3\pi/2","2\pi"]);
%xlim(ax,[klist_l(1),klist_l(end)]);

    %-------- return --------
    if nargout  == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    end
    if nargout  == 1
        varargout{1} = ax;
    end

end