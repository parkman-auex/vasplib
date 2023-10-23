function varargout = LocalChernplot(orb_list,WaveFunc,options)
arguments
    orb_list = [];
    WaveFunc = zeros(size(orb_list,1),1);
    options.ax =  handle([]);
    options.Rm = [];
    options.POSCAR = 'POSCAR';
    options.WaveMin = min(WaveFunc);
    options.WaveMax = max(WaveFunc);
%     options.WaveColor = 'r';
    options.WaveSize = 1;
    options.OrbColor = 'k';
    options.OrbSize = 1;
end

if isempty(options.Rm)
    try
        [Rm,~,~,~,~]=POSCAR_read(options.POSCAR);
    catch
        filename = input('please give a POSCAR file :');
        [Rm,~,~,~,~]=POSCAR_read(filename);
    end
else
    Rm = options.Rm;
end
if isempty(options.ax)
    Fig = Figs(1,1);
    ax = Fig.axes(1);
else
    if ishandle(options.ax)
        ax = options.ax;
    else
    end
end
[Nwave,Nlist] = size(WaveFunc);
Norb = length(orb_list);
if  Norb ~= Nwave
    error('Orbital list length is not equal to Local Chern list');
end
WFplot_list = zeros(Norb,4);
WFplot_list(:,1:3) = orb_list*Rm ;
WFplot_list(:,4) =WaveFunc;
%WFplot_list((WFplot_list(:,4)<options.WaveMin),4) = 0;
label_list = WFplot_list(:,4) >= options.WaveMin & WFplot_list(:,4) <= options.WaveMax;
% label_list2 = ~label_list;
% scatter3(ax,WFplot_list(label_list,1),WFplot_list(label_list,2),WFplot_list(label_list,3),(WFplot_list(label_list,4))*1000*options.WaveSize,'filled','MarkerFaceColor',options.WaveColor);

c=(WFplot_list(label_list,4));

scatter3(ax,WFplot_list(label_list,1),WFplot_list(label_list,2),WFplot_list(label_list,3),options.WaveSize*1000,c,'.');
hold(ax,'on');
grid(ax,'off');
plot3(ax,WFplot_list(:,1),WFplot_list(:,2),WFplot_list(:,3),'ko','MarkerSize',options.OrbSize,'MarkerFaceColor',options.OrbColor);

view(0,90);
axis(ax,'equal');
axis(ax,'off');

    if nargout  == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    end
    if nargout  == 1
        varargout{1} = ax;
    end
end
