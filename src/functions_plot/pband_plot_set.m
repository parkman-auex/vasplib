%% pband_plot_set: plot dos
%%
%% Usage:
%
% * [fig,ax] = pband_plot_set(EIGENCAR_list,WEIGHTCAR_list,Name_list,Ecut,titlestring)
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
% * Last updated : 2023/09/04
%
%% Source code :
%
function varargout = pband_plot_set(klist,EIGENCAR,WEIGHTCAR_cell,Name_list,Selected_projs,options,optionsplot)
arguments
    klist
    EIGENCAR
    WEIGHTCAR_cell
    Name_list
    Selected_projs double
    options.cmap
    optionsplot.density = 1;
    optionsplot.WEIGHTCAR_factor = 1;
    optionsplot.ax = handle([]);
    optionsplot.filled = true;
end
%--------  init  --------
Nbands = size(EIGENCAR,1);
klist_mat = repmat(klist,[Nbands,1]);
klist_vector = klist_mat(:);
EIGENCAR_vector = EIGENCAR(:);
%--------  narg  --------
if isempty(optionsplot.ax)
    Fig =  Figs(1,1);
    optionsplot.ax = Fig.axes(1);
else
    
end
%density = optionsplot.density;
%WEIGHTCAR_factor = optionsplot.WEIGHTCAR_factor;
%--------  chek  --------
if ~isa(WEIGHTCAR_cell,'cell')
    error('The data format of input WEIGHTCAR is wrong!');
end
nWEIGHTCAR = length(WEIGHTCAR_cell);
[nColor,~] = size(options.cmap);
if nColor < nWEIGHTCAR
    try
        cmap = options.cmap(nWEIGHTCAR);
    catch
        error('Color type is less than WEIGHTCAR');
    end
else
    cmap = options.cmap;
end
optionsplotcell = namedargs2cell(optionsplot);
%--------  plot  --------
for i = Selected_projs
    color = cmap(i,:);
    %disp(color);
    WEIGHTCAR_vector = WEIGHTCAR_cell{i}(:);
    displayname = Name_list(i,:);
    ax = pband_plot_one_bubble(klist_vector,EIGENCAR_vector,WEIGHTCAR_vector,color,displayname,optionsplotcell{:});
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