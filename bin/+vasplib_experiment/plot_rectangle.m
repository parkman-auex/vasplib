function [r_handle,f_handle] = plot_rectangle(pos,ax,color,name)
if nargin <3
    color = 'y';
end
if nargin <4
    name  = 'temp slice';
end
    r_handle = rectangle(ax,'Position',pos,'EdgeColor',[0.9,0.9,0.9],'LineWidth',0.2);
    Xlist = [pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)];
    Ylist = [pos(2),pos(2),pos(2)+pos(4),pos(2)+pos(4)];
    f_handle=fill(ax,Xlist, Ylist ,color,'facealpha',0.3);
end