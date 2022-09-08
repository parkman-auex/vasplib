function [parity]=Parity(filename,Occupiband)

% headline
%% Initialize variables.
startRow = 2;
%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%3s%3s%11s%11s%12s%[^\n\r]';
%% Open the text file.
fileID = fopen(filename,'r');
%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);
%% Allocate imported array to column variable names
bnd = double(dataArray{1});
nd = double(dataArray{2});

E = double(dataArray{4});
I = double(dataArray{5});
label = double(dataArray{6});
if nargin<2
    Efermi=textread('Efermi');
    geigval = double(dataArray{3})-Efermi;
    Occupiseq =find(geigval>0,1)-1;
else
    Occupiseq =find(bnd>Occupiband,1)-1;
end
parity.tot=sum(I(1:Occupiseq ));
 Minus  = 0;
 Plus = 0;
for i=1:Occupiseq
    if I(i ) >0

    Plus = Plus +I(i ) ;
    
    elseif I(i ) < 0
           Minus = Minus +I(i ) ;
    else
    end
end
if -Minus +Plus < Occupiband
    temp_num = (Occupiband+Minus-Plus)/2;
    Minus = Minus-temp_num;
    Plus = Plus+temp_num;
end
parity.minus = Minus;
parity.plus = Plus;
end