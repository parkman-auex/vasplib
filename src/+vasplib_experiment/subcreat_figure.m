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
function [fig,ax] = subcreat_figure(m,n,fontname,color,patersize,fig,ax)
%--------  init  --------
    FontSize = 24  ;
%--------  narg  --------
    if nargin <1
        m =2;
    end
    if nargin <2 
        n =1;
    end
    if nargin < 3
        fontname="Helvetica";
    end
    if nargin < 4
        color = 'white';
    end
    if nargin < 5
        patersize = [8 6];
    end
    if nargin < 6
        fig = figure('PaperType','a4letter','PaperSize',patersize,'Color',color);
    end
    if nargin < 7
        for i =1:m*n
            ax(i) = subplot(m,n,i,'Parent',fig ,'LineWidth',1.5,'FontSize',FontSize,'FontName',fontname); 
        end
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