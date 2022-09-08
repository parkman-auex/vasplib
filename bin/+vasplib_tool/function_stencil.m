%% pauli_matrice
%
% creat pauli_matrice: 
% $\sigma_0,\sigma_x,\sigma_y,\sigma_z$
%
% * Label:
%
%% Description of the Function:
%%
%% Usage: 
%
% * y-func (input1, input2);
% * y-func (input1, input2);
% * y-func (input1, input2);
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
function mat = pauli_matrice(dir)
%--------  init  --------
    I = 1i;
%--------  narg  --------

%--------  chek  --------
if nargin < 1
    sigma.o = eye(2);
    sigma.x = [ 0  1;...
                1  0];
    sigma.y = [ 0 -I;...
                I  0];
    sigma.z = [ 1  0;...
                0 -1];
    mat = sigma;
end
%--------  fbug  --------

%--------  juge  --------

%-------- return --------
end