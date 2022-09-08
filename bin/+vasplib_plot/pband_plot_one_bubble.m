%% pband_plot_one_bubble: plot pband
%%
%% Usage:
%
% * [fig,ax] = pband_plot_one_bubble(EIGENCAR_list,WEIGHTCAR_list,Name_list,Ecut,titlestring)
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
function varargout = pband_plot_one_bubble(klist,EIGENCAR,WEIGHTCAR,color,displayname,optionsplot)
arguments
    klist
    EIGENCAR
    WEIGHTCAR
    color
    displayname
    optionsplot.density = 1;
    optionsplot.WEIGHTCAR_factor = 1;
    optionsplot.ax = handle([]);
    optionsplot.filled = true;
end
%--------  init  --------
import vasplib_tool.*
%--------  init  --------
WEIGHTCAR = (WEIGHTCAR.*2.*10).^2;
WEIGHTCAR(WEIGHTCAR==0) = nan;
%Alpha_Data = mean(WEIGHTCAR_vector);
%--------  narg  --------
if isempty(optionsplot.ax)
    Fig =  Figs(1,1);
    ax = Fig.axes(1);
else
    ax = optionsplot.ax;
end

density = optionsplot.density;
WEIGHTCAR_factor = optionsplot.WEIGHTCAR_factor;
%--------  reshape data  --------
kn = length(klist);
kn_new = round(kn/density);
label_list = (1:kn_new)*density;
klist_new = klist(label_list);
EIGENCAR_vector_new = EIGENCAR(label_list);
WEIGHTCAR_vector_new = WEIGHTCAR(label_list);
WEIGHTCAR_vector_new = WEIGHTCAR_vector_new * WEIGHTCAR_factor;
%--------  main  --------
if optionsplot.filled
    scatter(ax,klist_new,EIGENCAR_vector_new,WEIGHTCAR_vector_new,'MarkerEdgeColor','none','MarkerFaceColor',color,'DisplayName',displayname);
else
    scatter(ax,klist_new,EIGENCAR_vector_new,WEIGHTCAR_vector_new,'MarkerEdgeColor',color,'MarkerFaceColor','none','DisplayName',displayname);
end
%--------  post  --------

%--------  save  --------
%-------- return --------
if nargout  == 2
    varargout{1} = ax.Parent;
    varargout{2} = ax;
end
if nargout  == 1
    varargout{1} = ax;
end
end