%% SOTI naïve model
%% 
% 
% 
% 2D Hcktpre model ta tb effect parkman 2022-01-12

% ------------------------    POSCAR    ------------------------
% TaTbHex
% 1.0
% 4.0000000000         0.0000000000         0.0000000000
% -2.0000000000         3.4641016151         0.0000000000
% 0.0000000000         0.0000000000         3.0000000000
% C Si
% 3  3
% Direct
% 0.300000012         0.000000000         0.000000000 C s 
% 0.000000000         0.300000012         0.000000000 C s 
% 0.699999988         0.699999988         0.000000000 C s 
% 0.699999988         0.000000000         0.000000000 Si s 
% 0.000000000         0.699999988         0.000000000 Si s 
% 0.300000012         0.300000012         0.000000000 Si s 
% ------------------------    ______    ------------------------
% %% Constuct TB model
% Type = 'list';%listmat
% HexTaTb = HR.from_POSCAR_SE('POSCAR','level_cut',2,'onsite',false,'per_dir',[1,1,0],'WAN_NUM',6,'r_max',3,'Type',Type);
% % 
% syms t_a t_b real;
% syms E1 E2 v w real;
% HexTaTb = HexTaTb.subs([E1 E2 t_a t_b]);
% HexTaTb_1 = subs(HexTaTb,[E1 E2 t_a t_b],[E1 -E1 v w]);
% HexTaTb_1 = HexTaTb_1<'POSCAR_sym';
%% armchair 2
Type = 'list';%listmat
syms v w real;
HexTaTb = HR.from_POSCAR_SE('POSCAR','level_cut',2,'onsite',false,'per_dir',[1,1,0],'WAN_NUM',4,'r_max',5,'Type',Type);
HexTaTb = subs(HexTaTb,[w v]);
HexTaTb = HexTaTb<'POSCAR';
%% prepare
Hcktpre = HexTaTb.HRforHckt();
%% Add Unit
L_0 = 1e-6;
MU = 100;
Cfactor = 100*1E-12;% 100 pF
OmegaFactor = (Cfactor*L_0*100)^(-1/2);
para_list = [1,3];

FigTest = Figs(1,2);
v= para_list(1,1);
w= para_list(1,2);
C_0 = 1;
% drawnow;
bandplot(HexTaTb.Subsall().EIGENCAR_gen('convention','I'),...
    [-5,5],'title',"t_b = "+num2str(w)+", t_a = "+num2str(v),'Color','b','ax',FigTest.axes(1));
%
Hcktpre_n = Hcktpre.Subsall();
%
EIGENCAR = Hcktpre_n.EIGENCAR_gen();
F0CAR = (EIGENCAR*L_0*Cfactor).^(-1/2)./(2*pi);
OMEGACUTpre = [0,2];
OMEGACUT = OMEGACUTpre*OmegaFactor;
%
bandplot(F0CAR,OMEGACUTpre*OmegaFactor,'ylabel','Frequency (Hz)','title',[...
    string(['C_A = ',num2str(MU*double(subs(Hcktpre.HcoeL(1)))),'pF']);...
    string(['C_v = ',num2str(MU*v),'pF','; C_w = ',num2str(MU*w),'pF']);...
    ],...
    'ax',FigTest.axes(2));