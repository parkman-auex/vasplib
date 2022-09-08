%% Parity TB test for TBG
%% 生成TBG TB
% addpath('..');
TBG = HR.from_POSCAR_SE('POSCAR',...
    'Type','mat',...
    'r_max',6,...
    'level_cut',-1,...
    'search_range',[1 1 0],...
    'Accuracy',1e-3,...
    'Rd',[sym('d0'),sym('a0')],...
    'E0',sym('E0'),...
    'chiral',false,...
    'onsite',false,...
    'deltarule',1);
TBG = TBG <'KPOINTS';
TBG.orbL
TBG.quantumL
TBG.symvar_list
%%  给出参数

% parm
E0 = 0.7632;
VppP_1  = 2.7;
VppS_5  = -0.48*1.7;
%Rnn
a0 = 1.42 ;
d0 = 3.35;  
delta=0.45255; % VppP

% 能带

TBG_n = TBG.Subsall();
EIGENCAR = TBG_n.EIGENCAR_gen();
bandplot(EIGENCAR ,[-3,3],'title',"TBG-TB-VppP_1="+string(VppP_1));
%% 调整TB轨道


%% 输出HR.dat

% TBG_n.Gen_hr('lda_hr.dat');
%% 自动给出TBmodel的输入文件
TBG_n.tbbox_in_gen();
%% check C2z
C2z = Oper.rotation(1/2,[0,0,1],false);
klist_TRIM = [0,0,0;0.5,0,0;0,0.5,0;0.5,0.5,0;];%0,0,0.5;0.5,0,0.5;0,0.5,0.5;0.5,0.5,0.5;];
[PARITYCAR,C2z] = C2z.Character_gen(TBG_n,klist_TRIM);
%% Wilsonloop
[BFCAR,~,klist_l] = TBG_n.WilsonLoop('knum_evol',101,'knum_int',51,'kstart',[0,-0.5,0.0],'kevolution',[1,0,0],'kintegral',[0,1,0]);
[fig,ax] = vasplib.WilsonLoopPlot(BFCAR,klist_l,'Color','r');
