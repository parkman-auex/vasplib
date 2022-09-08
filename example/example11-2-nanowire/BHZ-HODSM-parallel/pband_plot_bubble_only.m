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
% * Last updated : 2020/12/03
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code : 
%
function [fig,ax] = pband_plot_bubble_only(klist,EIGENCAR,WEIGHTCAR_cell,Name_list,cmap,fig,ax)
%--------  init  --------
    import vasplib_tool.*  
    Nbands = size(EIGENCAR,1);
    klist_mat = repmat(klist,[Nbands,1]);
    klist_vector = klist_mat(:);
    EIGENCAR_vector = EIGENCAR(:);
%--------  narg  --------
    if nargin < 6
         [fig,ax] = creat_figure();
    end
%--------  chek  --------
    if ~isa(WEIGHTCAR_cell,'cell')
        error('The data format of input WEIGHTCAR is wrong!');
    end
    nWEIGHTCAR = length(WEIGHTCAR_cell);
    [nColor,~] = size(cmap);
    if nColor < nWEIGHTCAR 
        error('Color type is less than WEIGHTCAR');
    end
%--------  plot  --------
for i = 1:nWEIGHTCAR
    %color = cmap(i,:);
    %disp(color);
    WEIGHTCAR_vector = WEIGHTCAR_cell{i}(:);
    displayname = Name_list(i,:);
    %--------  1 --------
    WEIGHTCAR_vector = (WEIGHTCAR_vector.*2.*10).^2;
    color_vector_init = normalize(WEIGHTCAR_vector,'range');
    color_vector_init(:,2:3) = ones(length(color_vector_init ),2);
    color_vector = hsv2rgb( color_vector_init);
    %WEIGHTCAR_vector(WEIGHTCAR_vector==0) = nan;
    %Alpha_Data = mean(WEIGHTCAR_vector);
    %--------  main  --------
    scatter(ax,klist_vector,EIGENCAR_vector,'CData',color_vector,'DisplayName',displayname);
    colormap(ax,cmap);
%     scatter(ax,klist_vector,EIGENCAR_vector,WEIGHTCAR_vector,'MarkerEdgeColor',color_vector,'MarkerFaceColor',color_vector,'DisplayName',displayname);
    %-------- 2--------
end
%--------  post  --------
%--------  save  --------
%-------- return --------
end