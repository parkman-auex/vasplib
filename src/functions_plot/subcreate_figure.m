%% subcreat_figure
%
% subcreat_figure for article print
%
% * Label:
%
%% Description of the Function:
%%
%% Usage:
%
% * [fig,ax] = subcreat_figure(fontname,color,patersize,fig,ax)
% * [fig,ax] = subcreat_figure(fontname,color,patersize)
% * [fig,ax] = subcreat_figure(fontname,color)
% * [fig,ax] = subcreat_figure(fontname)
% * [fig,ax] = subcreat_figure()
%
%% Input:
%
% # input1:
% # input2:
% # input3:
%
%% Output:
%
% # output1:
% # output2:
% # output3:
%
%% example:
%   commmad
%
%
%   result
%
%% Note:
%
%
%
%% Change log
%
% * Document Date: 2020/12/03
% * Creation Date: 2020/12/03
% * Last updated : 2020/12/03
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code :
%
function [fig,ax] = subcreate_figure(m,n,position,options)
arguments
    m double{mustBeInteger}=2;
    n double{mustBeInteger}=1;
    position = [];
    options.papersize = 1;
    options.fontname="Helvetica";
    options.color = 'white';
    options.fig = handle([]);
    options.ax  = handle([]);
    options.FontSize = 12;
end
if isempty(position)
    maxmn = max(m,n);
    if n == maxmn
        height  = 0.8* m/n;
        width = 0.6;
    end
    if m == maxmn
        height  = 0.8;
        width = 0.6*n/m;
    end

    position = [0.1 0.1 width height];
end
position([3,4]) = position([3,4])*options.papersize;
if isempty(options.fig )
    fig = figure('unit','normalized','position',position,'Color',options.color);
else
    fig =  options.fig ;
end

if isempty(options.fig )
    for i =1:m*n
        ax(i) = subplot(m,n,i,'Parent',fig ,'LineWidth',1.5,'FontSize',options.FontSize,'FontName',options.fontname);
    end
else
    ax = options.ax;
end

%--------  chek  --------

%--------  fbug  --------

%--------  main  --------
for i = 1:length(ax)
    box(ax(i),'on');
    hold(ax(i),'all');
end
%-------- return --------
end