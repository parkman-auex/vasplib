%% creat_figure
%
% creat_figure for article print
%
% * Label:
%
%% Description of the Function:
%%
%% Usage: 
%
% * [fig,ax] = creat_figure(fontname,color,patersize,fig,ax)
% * [fig,ax] = creat_figure(fontname,color,patersize)
% * [fig,ax] = creat_figure(fontname,color)
% * [fig,ax] = creat_figure(fontname)
% * [fig,ax] = creat_figure()
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
function [fig,ax] = creat_figure(fontname,color,patersize,fig,ax)
%--------  init  --------
    FontSize = 24  ;
%--------  narg  --------
    if nargin < 3
        patersize = [8 6];
    end
    if nargin < 2
        color = 'white';
    end
    if nargin < 1
        fontname="Helvetica";
    end
    if nargin < 4
        fig = figure('PaperType','a4letter','PaperSize',patersize,'Color',color);
    end
    if nargin < 5
        ax = axes('Parent',fig ,'LineWidth',1.5,'FontSize',FontSize,'FontName',fontname);  
    end
    
%--------  chek  --------

%--------  fbug  --------

%--------  main  --------
    box(ax,'on');
    hold(ax,'all');
%-------- return --------
end