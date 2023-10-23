%% dos_plot_set: plot dos 
%%
%% Usage: 
%
% * [fig,ax] = dos_plot_set(EIGENCAR_list,WEIGHTCAR_list,Name_list,Ecut,titlestring)
%
%% Input:
%  
% # input1:
% # input2:
% # input3:
%
%% Output:
%
% # fig:
% # ax:
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
%  Take advantage of the scope of application of the function.
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
function varargout = dos_plot_set(EIGENCAR_list,WEIGHTCAR_list,Name_list,Ecut,titlestring,options)
arguments
    EIGENCAR_list
    WEIGHTCAR_list
    Name_list
    Ecut
    titlestring
    options.cmap = @jet;
end
%--------  init  --------
import vasplib_tool.*
%--------  init  --------
[~,ax] = create_figure('Position',[0.2,0.2,0.2,0.6]);
%--------  narg  --------

%--------  chek  --------
%sColorMap;
cmap = options.cmap(length(Name_list));
%--------  plot  --------
[Nname,~]=size(Name_list);
for i =1:Nname
    plot(ax,WEIGHTCAR_list(:,i),EIGENCAR_list(:,i),'Color',cmap(i,:),'LineWidth',2.0,'DisplayName',Name_list(i,:));
end
%--------  post  --------
ax.YLim = Ecut;
ylabel(' E-E_f (eV) ');
legend();
title(titlestring);
%--------  save  --------
%[fig,ax] = save_figure(fig,titlestring+".eps",ax);
%-------- return --------
if nargout  == 2
    varargout{1} = ax.Parent;
    varargout{2} = ax;
end
if nargout  == 1
    varargout{1} = ax;
end
end