function [DOSCAR,X,Y] = DOSCAR_read(filename,YXlength,mode,options)
%IMPORTFILE 将文本文件中的数值数据作为矩阵导入。
%   [DOSCAR,X,Y]  = DOSCAR_read(FILENAME)
%   读取文本文件 FILENAME 中默认选定范围的数据。
%
%   [DOSCAR,X,Y]  = DOSCAR_read(FILENAME,'surf',STARTROW, ENDROW)
%   读取文本文件 FILENAME 的 STARTROW 行到 ENDROW 行中的数据。
%
% Example:
%   dos1 = DOSCAR_read('dos.dat_r');
%
%    另请参阅 TEXTSCAN。

% 由 MATLAB 自动生成于 2021/12/07 21:19:21
%% 初始化变量。
arguments
    filename = 'dos.dat_r';
    YXlength = [101 1000];
    mode = 'surf';
    options.startRow = 1
    options.endRow = inf
end
% startRow = options.startRow ;
% endRow = options.endRow ;
switch mode
    case 'surf'
        DOSCAR_raw = importdata(filename);
    case 'arc'
        DOSCAR_raw = importdata(filename);
        DOSCAR_raw = DOSCAR_raw.data;
end
if isa(DOSCAR_raw,'struct')
    %DOSCAR_raw = importdata(filename);
    DOSCAR_raw = DOSCAR_raw.data;
end
%% 创建输出变量
X = reshape(DOSCAR_raw(:,1),YXlength).';
Y = reshape(DOSCAR_raw(:,2),YXlength).';
DOSCAR = reshape(DOSCAR_raw(:,3),YXlength).';