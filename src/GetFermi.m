%% Get fermi
%
%
% * Label: Tools function
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
function Efermi = GetFermi(mode,filename)
if nargin < 2
    filename = 'DOSCAR';
end
    if strcmp(mode,'vasp')
        % read fermi
        
        delimiter = ' ';
        startRow = 6;
        endRow = 6;
        formatSpec = '%*s%*s%*s%f%*s%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        Efermi = [dataArray{1:end-1}];
        clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;
    elseif strcmp(mode,'qe')
        [Efermi_string,~]=grep('scf.out','Fermi');
        Efermi = double(regexp(Efermi_string,'\d*\.?\d*','match'));
        disp(Efermi);
    else
        Efermi = 0 ;
    end

end