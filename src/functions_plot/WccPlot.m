function varargout = WccPlot(WccCAR,klist_l,options)
arguments
    WccCAR;
    klist_l;
    options.ycut = [0,1];
    options.shift = 0;
    options.Color = [rand rand rand];
    options.LineSpec = 'o';
    options.title = '';
    options.xlabel='k_{envolution}';
    options.ylabel=["Wcc"];
    options.ax = handle([]);
end
if options.ycut(1) >= 0
    WccCAR(WccCAR(:,:)<0) = 1+WccCAR(WccCAR(:,:)<0);
end
if isempty(options.ax)
    Fig = create_figure('Position',[0.2,0.2,0.6,0.6]);
    ax = Fig.axes(1);
else
    ax = options.ax ;
end
if length(klist_l)==1
    WccCAR = [WccCAR,WccCAR];
    klist_l = [0,klist_l];
end
ax = bandplot(WccCAR+options.shift,options.ycut+options.shift,klist_l,[],[],...
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
set(ax,'ytick',-1:0.25:1,...
    'yticklabel',["-1","-0.75","-0.5","-0.25",...
    "0","0.25","0.5","0.75","1"]);
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