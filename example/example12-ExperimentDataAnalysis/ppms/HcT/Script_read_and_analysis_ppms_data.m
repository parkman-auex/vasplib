%%
% A fast script for PPMS hc data import  
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
ISTART = 3 ; % 0,1,2,3 for different start level
Sample_name = 'NaYbO2' ; % Give a fancy name for your sample
Sample_mass = 0.0049   ;%g
Sample_RFM  = 223.89      ;
Sample_num  = Sample_mass/Sample_RFM;% n mol
Data_filename = '20201219-Sample 4.9mg.Dat' ;
Data_density = 100     ; % for plot, 100 ~ 150 recommand
%-------------------------------------------------------------------------------
%% CONTROL
fit_Shutty = 1 ; 
Trange_list = [59,81;
               59,81;
               59,81]; %K for CW
Trange_plot = [ 0,20;... 
                0,20;...
                0,20;...
                ]; %K for plot
C_sch_input.n = Sample_num ;              % mol 
%C_sch_input.R = 8.314472  ;               % J/(K * mol)
                        ND0 = 1;ND1 = 1;
C_sch_input.DegenRate = ND0/ND1;
C_sch_input.Delta = 0.99*0.53 ;           % meV
% for fit
T_range_extend = [2,20;20,60;65,300];

Low_value   =   [ 1e-5  1e-8 1e-16 30];
% Start_point =   ;
High_value  =   [ 1e-2 1e-3 1e-3 100];
%-------------------------------------------------------------------------------
%% ISTART
%-------------------------------------------------------------------------------
%********************************************************************************
%% From different start level
        
switch ISTART
    case 0
        PPMS_HC_DATA= importPpmsHcData(Data_filename);
        [HcDesignData,H_total,T_total,Hc_total,time_total] = PpmsHcData_gen(PPMS_HC_DATA,Sample_num );
        disp(PPMS_HC_DATA(1,:));
% Take Hc H T data
% in Low temperature Hc measurement, we focus on the data of Hc, H, T, 
% explicitly as (\muJ/K), MagField(Oe) and Temperature(K)
        disp(HcDesignData(1,:));
        [Hc_cell,T_cell,H_cell,H_list] = Hc_T_split(HcDesignData);
        save('PPMS_HC_DATA.mat','PPMS_HC_DATA','HcDesignData',...
            'H_total','T_total','Hc_total','time_total',...
            'Hc_cell','T_cell','H_cell','H_list');   
    case 1
        load PPMS_HC_DATA.mat;
    case 2
        load PPMS_HC_DATA.mat;
        disp('we will not choose the T_range for DebyeT3');
    case 3
        load PPMS_HC_DATA.mat;
        disp('we will fit Shutty');
end

%% Data plot
% ------------------------------- t3 fit ------------------------------
% we offer two mode : with \chi_0 and without \chi_0
%
if ISTART  == 1

    nSituation= length(H_list);
    for i = 1:nSituation
        % Data
        Hc = Hc_cell{i};
        T = T_cell{i};
        Name = "H = "+string(H_list(i))+" Ost";
        % fit       
        [FittingData(i),Hc_brandnew{i},fig(i),ax(i)]=Debye_T3(Hc,T,Data_density,Name);
    end
    % ------------------------------- HcH plot ------------------------------
    %
    [fig,ax] = HcT_figure_plot(Hc_brandnew,T_cell,H_list,Data_density,Trange_plot);

elseif ISTART  == 2

    nSituation= length(H_list);
    for i = 1:nSituation
        % Data
        Hc = Hc_cell{i};
        T = T_cell{i};
        Name = "H = "+string(H_list(i))+" Ost";
        T_range = Trange_list(i,:);
        % fit
        [FittingData(i),Hc_brandnew{i},fig(i),ax(i)]=Debye_T3(Hc,T,Data_density,Name,T_range);
    end
    % ------------------------------- HcH plot ------------------------------
    %
    [fig,ax] = HcT_figure_plot(Hc_brandnew,T_cell,H_list,Data_density,Trange_plot);
elseif ISTART  == 3
    nSituation= fit_Shutty;
    % Data
    Hc = Hc_cell{nSituation};
    T = T_cell{nSituation};
    Name = "H = "+string(H_list(nSituation))+" Ost";
    T_range = Trange_list(nSituation,:);
    
    % fit
    [FittingData,Hc_brandnew,fig,ax]=LT_specific_heat(Hc,T,Data_density,Name,C_sch_input,...
        T_range_extend,Low_value,High_value);

end

%% ask if save the workspace

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
