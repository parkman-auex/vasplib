%% bandplot
%
% plotband
%
%
% * Label: plot
%
%% Description of the Function:
%%
%% Usage:
%
% * [fig,ax]=bandplot(EIGENCAR,Ecut,titlestring,color,klist_l,kpoints_l,kpoints_name,fontname,fig,ax)
% * [fig,ax]=bandplot(EIGENCAR,Ecut,titlestring,color,klist_l,kpoints_l,kpoints_name,fontname)
% * [fig,ax]=bandplot(EIGENCAR,Ecut,titlestring,color,klist_l,kpoints_l,kpoints_name)
% * [fig,ax]=bandplot(EIGENCAR,Ecut,titlestring,color)
% * [fig,ax]=bandplot(EIGENCAR,Ecut,titlestring)
% * [fig,ax]=bandplot(EIGENCAR,Ecut)
% * [fig,ax]=bandplot(EIGENCAR)
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
%
%% example:
% for vasp
%   bandplot();
% for a eigencar
%   bandplot(EIGENCAR)
% use Ecut
%   bandplot(EIGENCAR,Ecut)
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
function varargout=bandplot(varargin)
[varargout{1:nargout}] = vasplib_plot.bandplot(varargin{:});
end


