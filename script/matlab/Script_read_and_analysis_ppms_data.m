%%
% A fast script for PPMS data import  
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Let's begin
% addpath('~/Documents/vasplib');
%% As a envionment start, you should load envionment
import vasplib_experiment.*
%-------------------------------------------------------------------------------
%% SYMTEM
ISTART = 3; % 1,2,3 for different start level
Sample_name = 'NaTmO2' ; % Give a fancy name for your sample
Sample_mass = 0.1442   ;%g
Sample_RFM  = 142      ;
Sample_num  = Sample_mass/Sample_RFM;% n 
Data_filename = '20201202-NaTmO2-OP.DAT' ;
Data_density = 100     ;% for plot, 30 ~ 50 recommand
%-------------------------------------------------------------------------------
%% CONTROL
MH_seq_list = [5,6,7];
MH_T_list = [2,3,5]; % length equals to MH_seq_list
ZFC_FC_seq_list = [1,2;...
                   3,4]; 
ZFC_FC_name_list = ["100 Oe";"1T "];
Trange_CW = [100,300]; %K for CW
Trange_plot = [2,300]; %K for plot
Hrange = [0,8]*1000; %T
%-------------------------------------------------------------------------------
%% ISTART
%
% * for ISTART = 1
Block_list = [1,3369 ;3833 , 7541;  7870 , 11698; 11698 , 15621;15734 , 16132;15734 , 16132;15734 , 16132];
% * for ISTART = 2
%
% If your are familar to the code and data, we recommand using ISTART = 2
Block_struct(1).block = [1,3369];
    Block_struct(1).type = 2;
        Block_struct(1).name = '100_Oe';
%
% * for ISTART = 3
%
% make sure you have 'PPMS_VSM_DATA.mat'
%-------------------------------------------------------------------------------
%********************************************************************************
%% From different start level
        PPMS_VSM_DATA= importPpmsVsmData(Data_filename);
        [MagMeasureData,M_total,H_total,T_total,time_total] = PpmsMagData_gen(PPMS_VSM_DATA);
switch ISTART
    case 0
        disp(PPMS_VSM_DATA(1,:));
% Take M H T data
% in Low temperature magnetic measurement, we focus on the data of M, H, T, 
% explicitly as Magnetic Momentum(emu.), MagField(Oe) and Temperature(K)
        disp(MagMeasureData(1,:));
        [Block_struct,~,~,~,~] = Changelabel_by_HTpeaks(time_total,T_total,H_total);
        VsmMagDataStruct = VsmMagDataStruct_gen(Block_struct,MagMeasureData,Sample_name,PPMS_VSM_DATA);
        save('PPMS_VSM_DATA.mat');
    case 1
        [Block_struct,~,~,~,~] = Changelabel_by_HTpeaks(time_total,T_total,H_total,Block_list);
        VsmMagDataStruct = VsmMagDataStruct_gen(Block_struct,MagMeasureData,Sample_name,PPMS_VSM_DATA);
        save('PPMS_VSM_DATA.mat');
    case 2
        VsmMagDataStruct = VsmMagDataStruct_gen(Block_struct,MagMeasureData,Sample_name,PPMS_VSM_DATA);
        save('PPMS_VSM_DATA.mat');
    case 3
        load('PPMS_VSM_DATA.mat');
end
%% Data plot
% ------------------------------- Cuire Weiss fit ------------------------------
% we offer two mode : with \chi_0 and without \chi_0
%
[nSituation,~] = size(ZFC_FC_seq_list);
for i = 1:nSituation
    % ZFC
    [M_ZFC{i},H_ZFC{i},T_ZFC{i}] = ...
        MHT_from_table(VsmMagDataStruct(ZFC_FC_seq_list(i,1)).datatable);
    % FC
    [M_FC{i},H_FC{i},T_FC{i}] = ...
        MHT_from_table(VsmMagDataStruct(ZFC_FC_seq_list(i,2)).datatable); 
    % ZFC
    title = strcat(ZFC_FC_name_list(i,:),'ZFC');
    [FittingData1_ZFC_temp,fig1_CW(i),ax1_CW(i)]=...
        Curie_Weiss(M_ZFC{i},H_ZFC{i},T_ZFC{i},Sample_num,Trange_CW,Data_density,title);
    FittingData1_ZFC{i} = FittingData1_ZFC_temp;
    % FC
    title = strcat(ZFC_FC_name_list(i,:),'FC');
    [FittingData2_FC_temp,fig2_CW(i),ax2_CW(i)]=...
        Curie_Weiss(M_FC{i},H_FC{i},T_FC{i},Sample_num,Trange_CW,Data_density,title); 
    FittingData1_FC{i} = FittingData2_FC_temp;
    % 
    title = ZFC_FC_name_list(i,:);
    [fig_chichim1T(i),ax_chichim1T(i)] = ...
        X_and_Xm1_figure_plot(M_ZFC{i},H_ZFC{i},T_ZFC{i},...
                               M_FC{i}, H_FC{i}, T_FC{i},Sample_num,Trange_plot,Data_density,title);
end
% ------------------------------- MH plot ------------------------------    
%
[M_cell,H_cell,T_cell] = MHT_cell_from_table_struct(VsmMagDataStruct,MH_seq_list);
[fig_MH,ax_MH] = MH_figure_plot(M_cell,H_cell,MH_T_list,Data_density);
%% ask if save the workspace

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
