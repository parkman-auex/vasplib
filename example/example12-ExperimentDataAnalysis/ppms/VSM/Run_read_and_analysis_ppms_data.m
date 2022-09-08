%% example12 ExperimentDataAnalysis
% ****************************************
%             # ppms ssm data read and analysis
% ****************************************
%
%% Label:
% * ppms
% * mag
% * data cleam
%
%% Description of the Script:
%
% make a ppms vsm data analysis
%
%% Usage: 
%
% * Modify or Understand this code
%
%% Input:
%  
% # Files: a PPMS_VSM.DAT
%
%% Output:
%
% # Typical serveral data and plot
%
%% Note: 
%
%  Take advantage of the scope of application of the function.
%
%% Change log
%
% * Document Date: 2020/12/06
% * Creation Date: 2020/12/06
% * Last updated : 2020/12/08
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
% Let's begin
%
%% make sure you have the PPMS_VSM gen data in this directory
ls

%% As a envionment start, you should load envionment
import vasplib_experiment.*

%% for a default print,just type 
%   ppms()
% or 
%   ppms
%ppms;

%% for a user-defined ppms data analysis
%
% we shoulf first load the PPMS data, 
%
%   PPMS_VSM_DATA= importPpmsVsmData('PpmsVsm.DAT');
%
% which is a highly useful Tabel stucture
% with the format
%
% <html>
% <style type="text/css">
% .tg  {border-collapse:collapse;border-spacing:0;}
% .tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
%   overflow:hidden;padding:10px 5px;word-break:normal;}
% .tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
%   font-weight:normal;overflow:hidden;padding:10px 5px;word-break:normal;}
% .tg .tg-0pky{border-color:inherit;text-align:left;vertical-align:top}
% </style>
% <table class="tg">
% <thead>
%   <tr>
%     <th class="tg-0pky">Time_Stamp_sec</th>
%     <th class="tg-0pky">Temperature_K</th>
%     <th class="tg-0pky">Magnetic_Field_Oe</th>
%     <th class="tg-0pky">Moment_emu</th>
%     <th class="tg-0pky">M_Std_Err_emu_</th>
%     <th class="tg-0pky">...</th>
%   </tr>
% </thead>
% <tbody>
%   <tr>
%     <td class="tg-0pky">.</td>
%     <td class="tg-0pky">230</td>
%     <td class="tg-0pky">99.83</td>
%     <td class="tg-0pky">0</td>
%     <td class="tg-0pky">0</td>
%     <td class="tg-0pky"></td>
%   </tr>
%   <tr>
%     <td class="tg-0pky">...</td>
%     <td class="tg-0pky">...</td>
%     <td class="tg-0pky">...</td>
%     <td class="tg-0pky">...</td>
%     <td class="tg-0pky">...</td>
%     <td class="tg-0pky">...</td>
%   </tr>
% </tbody>
% </table>
% </html>
%
% in this example, we take NaTmO2 as a trial

PPMS_VSM_DATA= importPpmsVsmData('20201202-NaTmO2-OP.DAT');
disp(PPMS_VSM_DATA(1,:));

%% Take M H T data
% in Low temperature magnetic measurement, we focus on the data of M, H, T, 
% explicitly as Magnetic Momentum(emu.), MagField(Oe) and Temperature(K)
 [MagDesignData,M,H,T,time] = PpmsMagData_gen(PPMS_VSM_DATA);
 disp(MagDesignData(1,:));
 
%% Data clean
% we always get a redundacy data. Clean it first
[Block_struct,~,~,~,~] = Changelabel_by_HTpeaks(time,T,H);
% gen different data in your directory for next use
VsmMagDataStruct = VsmMagDataStruct_gen(Block_struct,MagDesignData);
%% 
