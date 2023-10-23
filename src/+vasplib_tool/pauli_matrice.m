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
% * sigma = pauli_matrice();
% * sigma_i = pauli_matrice(dir_i);
%
%% Input:
%  
% # dir : 0;1;2;3;o;x;y;z
%
%% Output:
%
% # output1: sigma_mat(sigma_mat.o,sigma_mat.x,sigma_mat.y,sigma_mat.z)
% # output2: sigma_i
%
%% example:
%   tau = pauli_matrice();
%   tau.o;
%   tau.x;
%   tau.y;
%   tau.z;
%   
%% Note: 
%
%  Nothing.
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
if nargin < 1
    sigma.o = eye(2);
    sigma.x = [ 0  1;...
                1  0];
    sigma.y = [ 0 -I;...
                I  0];
    sigma.z = [ 1  0;...
                0 -1];
    mat = sigma;
    return;
end
%--------  chek  --------

%--------  fbug  --------

%--------  juge  --------
    %input_char_data = class(dir);
    input_char = char(string(dir));
    switch input_char
        case {'0','o'}
            mat = eye(2);
        case {'1','x'}
            mat =  [ 0  1;...
                    1  0];
        case {'2','y'}
            mat = [ 0 -I;...
                    I  0];
        case {'3','z'}
            mat = [ 1  0;...
                    0 -1];
    end    
%-------- return --------
end