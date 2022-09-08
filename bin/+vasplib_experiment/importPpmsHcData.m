function PPMS_HC_DATA= importPpmsHcData(filename, startRow, endRow)
%IMPORTFILE 从文本文件中导入数据
%  SAMPLE4 = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  以表形式返回数据。
%
%  SAMPLE4 = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  Sample4 = importfile("D:\zcdxsh\比热数据\GJ\20201219-Sample 4.9mg.Dat", [16, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2020-12-22 12:22:07 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
     startRow = 16;
end
if nargin < 3
     endRow = inf;
end
dataLines = [ startRow, endRow];
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 30);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = ",";

% 指定列名称和类型
opts.VariableNames = ["Time_Stamp_sec", "Comment", "System_Status_Code", "Puck_Temp_K", "System_Temp_K", "Magnetic_Field_Oe", "Pressure_Torr", "Temperature_K", "Temp_Rise_K", "SampleHC_muJ_d_K", "Samp_HC_Err_K", "Addenda_HC_K", "Addenda_HC_Err_K", "Total_HC_K", "Total_HC_Err_K", "Fit_Deviation_Chi_Square", "Time_Const_Tau1_seconds", "Time_Const_Tau2_seconds", "Sample_Coupling_Percent", "Debye_Temp_K", "Debye_Temp_Err_K", "Cal_Correction_Factor", "Therm_Resist_Ohms", "Htr_Resist_Ohms", "Puck_Resist_Ohms", "Wire_Cond_W_d_K", "Meas_Time_seconds", "Temp_Squared_K_2", "Samp_HC_Temp_KK", "Addenda_Offset_HC_K"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "double", "string", "string", "string", "double", "double", "double", "double", "double"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, ["Comment", "Debye_Temp_K", "Debye_Temp_Err_K", "Therm_Resist_Ohms", "Htr_Resist_Ohms", "Puck_Resist_Ohms"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Comment", "Debye_Temp_K", "Debye_Temp_Err_K", "Therm_Resist_Ohms", "Htr_Resist_Ohms", "Puck_Resist_Ohms"], "EmptyFieldRule", "auto");

% 导入数据
PPMS_HC_DATA = readtable(filename, opts);

end