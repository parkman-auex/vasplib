%% set_reference
%
% set_reference for article print
%
% * Label:
%
%% Description of the Function:
%%
%% Usage: 
%
% * [fig,ax] = set_reference(fig,figname,ax)
% * [fig,ax] = set_reference(fig,figname)
% * [fig,ax] = set_reference(fig)
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
function varargout = set_reference(kpoints_l,kpoints_name,X_cut,Y_cut,options)
%--------  init  --------
arguments 
    kpoints_l = [];
    kpoints_name = [];
    X_cut = [];
    Y_cut = [];
    options.mode  = 'band';
    options.ax  = handle([]);
    options.FontName = 'Helvetica';
    options.FontSize = 24;
    options.xlabel='';
    options.ylabel='E(eV)';
end
%--------  chek  --------
FontName = options.FontName;
%--------  fbug  --------
if isempty(options.ax)
    ax = gca;
else
    ax = options.ax;
end
%--------  main  --------
if strcmp(options.mode,'band')
% set label
try
    set(ax,'XLim',X_cut);
catch
    warning('check if load KPOINTS./Or use bandplot wrongly')
end
    set(ax,'YLim',Y_cut);
    %set(ax,'Xlabel',options.xlabel)  ;
    set(ax,'FontName',options.FontName,'FontSize',options.FontSize,'LineWidth',1,...
    'XTick',kpoints_l,'XTickLabel',...
    kpoints_name);
    ylabel(ax,options.ylabel,'FontName',FontName);
    xlabel(ax,options.xlabel,'FontName',FontName);
% reference line    
    %line(X_cut,[0 0],'LineWidth',0.1,'Color',[0. 0. 0.],'DisplayName','fermi-level')
    for i=1:(length(kpoints_name)-2)                                           
        X=[kpoints_l(i+1) kpoints_l(i+1)];
        line(ax,X,Y_cut,'LineWidth',0.1,'Color',[0 0. 0.],'DisplayName','K-path');
    end
end
hold(ax,'on');
axis(ax,'square');
%-------- return --------
if nargout  == 2
    varargout{1} = ax.Parent;
    varargout{2} = ax;
end
if nargout  == 1
    varargout{1} = ax;
end
%set(fig,'InvertHardCopy','off');
end