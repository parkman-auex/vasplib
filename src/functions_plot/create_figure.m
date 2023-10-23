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
function varargout = create_figure(M,N,propArgs,FigureOption)
    arguments
        M = 1;
        N = 1;
        propArgs.test = [];
%         propArgs.?Figs;
        % https://ww2.mathworks.cn/help/matlab/matlab_prog/function-argument-validation-1.html#mw_1b62b6d6-a445-4c55-a9b9-9c70becfdbe6      
        FigureOption.Position = [];
    end
    % 创建 figure
    propArgsCell = namedargs2cell(propArgs);
    Fig = Figs(M,N);
    %--------  init  --------
    if isempty(FigureOption.Position)
    else
        set(Fig.fig,'Position',FigureOption.Position);
    end
    %-------- return --------
    if nargout  == 2
        varargout{1} = Fig.fig;
        varargout{2} = Fig.axes;
    end
    if nargout  == 1
        varargout{1} = Fig;
    end
end