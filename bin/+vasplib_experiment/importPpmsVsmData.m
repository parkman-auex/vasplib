function PPMS_VSM_DATA = importPpmsVsmData(filename, startRow, endRow)
%IMPORTFILE 将文本文件中的数值数据作为矩阵导入。
%   PPMS_VSM_DATA = importPpmsVsmData(filename) 
%   读取文本文件 FILENAME 中默认选定范围的数据。
%   PPMS_VSM_DATA = importPpmsVsmData(filename, startRow) 
%   读取文本文件 FILENAME 的 STARTROW   行到 最后一 行中的数据。
%   PPMS_VSM_DATA = importPpmsVsmData(filename, startRow, endRow) 
%   读取文本文件 FILENAME 的 STARTROW   行到 ENDROW 行中的数据。
%
% Example:
%   NaTmO2_DATA = importPpmsVsmData('20201202-NaTmO2-OP.DAT', 36, 21759);
%
%   
import vasplib_experiment.*;
%% 初始化变量。
delimiter = ',';
if nargin<=2
    startRow = 36;
    endRow = inf;
end
formatSpec = '%*q%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
%% 打开文本文件。
fileID = fopen(filename,'r');
%% 根据格式读取数据列。
% 该调用基于生成此代码所用的文件的结构。如果其他文件出现错误，请尝试通过导入工具重新生成代码。
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% 关闭文本文件。
fclose(fileID);

%% 对无法导入的数据进行的后处理。
% 在导入过程中未应用无法导入的数据的规则，因此不包括后处理代码。要生成适用于无法导入的数据的代码，请在文件中选择无法导入的元胞，然后重新生成脚本。

%% 创建输出变量
PPMS_VSM_DATA  = table(dataArray{1:end-1}, 'VariableNames', {'Time_Stamp_sec','Temperature_K','Magnetic_Field_Oe','Moment_emu','M_Std_Err_emu_','Transport_Action','Averaging_Time_sec','Frequency_Hz','Peak_Amplitude_mm','Center_Position_mm','Coil_Signal_mV','Coil_Signa_mV','Range_mV','M_Quad_Signal_emu','M_Raw_emu','M_Raw_2emu','Min_Temperature_K_','Max_Temperature_K_','Min_Field_Oe_','Max_Field_Oe_','Mass_grams_','MotorLag_deg_','Pressure_Torr_','VSMStatus_code','Motor_Status_code','Measure_StaMeasureStatus_code_tus_code','Measure_Count','System_Temp__K','Temp_Status_code','FieldStatus_code','ChamberStatus_code','MotorCurrent_amps_','MotorHeatsink_Temp_C'});
end
