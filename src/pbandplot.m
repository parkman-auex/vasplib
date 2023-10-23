%% pbandplot
%
% plotpband 
% 
%
% * Label: plot
%
%% Description of the Function:
%%
%% Usage: 
%
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring,cmap,klist_l,kpoints_l,kpoints_name,fontname,fig,ax)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring,cmap,klist_l,kpoints_l,kpoints_name,fontname)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring,cmap,klist_l,kpoints_l,kpoints_name)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring,cmap)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct)
% * [fig,ax]=pbandplot()
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
% * Document Date: 2020/12/04
% * Creation Date: 2020/12/04
% * Last updated : 2020/12/04
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code : 
%
function varargout=pbandplot(varargin)
[varargout{1:nargout}] = vasplib_plot.pbandplot(varargin{:});
end