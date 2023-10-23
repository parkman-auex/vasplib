%% save_figure
%
% save_figure for article print
%
% * Label:
%
%% Description of the Function:
%%
%% Usage: 
%
% * [fig,ax] = save_figure(fig,figname,ax)
% * [fig,ax] = save_figure(fig,figname)
% * [fig,ax] = save_figure(fig)
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
function [fig,ax] = save_figure(fig,figname,ax)
     mkdir('plot_results');
%--------  init  --------
    figname=strrep(figname,' ','');
    figname=strrep(figname,'*','_multi_');  
    figname=strrep(figname,'-','_m_'); 
    figname=strrep(figname,'+','_p_');
    figname=strrep(figname,'/','_divid_');
    figname=strrep(figname,'~','_power_');    
%--------  narg  --------

    if nargin < 2
        ax = gca(fig);
    end
    if nargin <1
        % Create title
        dirname=pwd;
        dirname=strsplit(dirname,'/');
        figname=dirname(length(dirname));
        figname = strcat(figname,'.eps');
    end
    
%--------  chek  --------

%--------  fbug  --------

%--------  main  --------
    set(ax,'color','none');
    set(fig,'color','none');
    set(fig,'InvertHardCopy','off');
    %saveas(figure_k,titlename);
    % disp(class(fig));
    % disp(figname);
    figname = "plot_results/"+figname;
    print(fig,figname,'-depsc','-tiff','-cmyk');
%-------- return --------
    set(ax,'color','white');
    set(fig,'color','white');
    %set(fig,'InvertHardCopy','off');
end