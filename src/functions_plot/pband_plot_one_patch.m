%% pband_plot_one_patch: plot dos 
%%
%% Usage: 
%
% * [fig,ax] = pband_plot_one(EIGENCAR_list,WEIGHTCAR_list,Name_list,Ecut,titlestring)
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
function varargout = pband_plot_one_patch(klist,EIGENCAR,WEIGHTCAR,options)
arguments
    klist
    EIGENCAR
    WEIGHTCAR
    options.cmap
    options.ax = handle([]);
    options.LineWidth = 2;
end

%--------  init  --------
    Nbands=size(EIGENCAR,1);
    cmap = options.cmap;
%--------  narg  --------
    if isempty(options.ax)
        ax= create_figure();
    else
        ax = options.ax;
    end
%--------  main  --------
    for Ei=1:Nbands
        patch(ax,[klist NaN],[EIGENCAR(Ei,:) NaN],[WEIGHTCAR(Ei,:) NaN],'Marker','none','EdgeColor','interp','MarkerFaceColor','flat','LineWidth',options.LineWidth,'DisplayName',num2str(Ei));
        %hold on
    end
%--------  post  --------
    colormap(ax,cmap);
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