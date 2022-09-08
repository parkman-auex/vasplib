%% 创建一个KramerSSHTB模型

KramerSSH = HR.from_POSCAR_SE('POSCAR',...
'Type','list',...
'onsite',false,...
'Chiral',true,...
'search_range' ,[1 0 0],...
'r_max',5,...
'level_cut',5);
KramerSSH.list;
KramerSSH = KramerSSH <'POSCAR';
%%
syms C_inter v w Minus real;
KramerSSH = KramerSSH.subs([v w C_inter]);
Gauge = 1;
switch Gauge
    case 1
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 4 & KramerSSH.vectorL(:,5) == 1) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 1 & KramerSSH.vectorL(:,5) == 4) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 5 & KramerSSH.vectorL(:,5) == 8) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 8 & KramerSSH.vectorL(:,5) == 5) = Minus*C_inter;
    case 2
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 6 & KramerSSH.vectorL(:,5) == 7) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 7 & KramerSSH.vectorL(:,5) == 6) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 5 & KramerSSH.vectorL(:,5) == 8) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 8 & KramerSSH.vectorL(:,5) == 5) = Minus*C_inter;
    case 3
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 2 & KramerSSH.vectorL(:,5) == 3) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 3 & KramerSSH.vectorL(:,5) == 2) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 5 & KramerSSH.vectorL(:,5) == 8) = Minus*C_inter;
        KramerSSH.HcoeL(KramerSSH.vectorL(:,4) == 8 & KramerSSH.vectorL(:,5) == 5) = Minus*C_inter;
end

%QTI = QTI.autohermi();

HcktTB = KramerSSH;
%%

%% Para
Minus = -1;
C_inter = 1;
v = 1;
w = 2;
C_0 = 1;
% Zhao prb
% t1 = 1e-4; % transverse bond
% 
% t2 = 1*1e-4;
% t3 = 4*1e-4;
% 
% m1 = -1*1e-4;
% m2 = -4e-4;





% Type = 'list';%listmat
% syms v w real;
% HexTaTb = HR.from_POSCAR_SE('POSCAR','level_cut',2,'onsite',false,'per_dir',[1,1,0],'WAN_NUM',4,'r_max',5,'Type',Type);
% HexTaTb = subs(HexTaTb,[w v]);
% HexTaTb = HexTaTb<'POSCAR';
%% prepare
Hcktpre = HcktTB.HRforHckt();
%% Add Unit
L_0 = 1e-6;
MU = 100;
Cfactor = 100*1E-12;% 100 pF
OmegaFactor = (Cfactor*L_0*100)^(-1/2);

FigTest = Figs(1,2);


% drawnow;
bandplot(HcktTB.Subsall().EIGENCAR_gen('convention','I'),...
    [-5,5],'Color','b','ax',FigTest.axes(1));
%
Hcktpre_n = Hcktpre.Subsall();
%
EIGENCAR = Hcktpre_n.EIGENCAR_gen();
F0CAR = (EIGENCAR*L_0*Cfactor).^(-1/2)./(2*pi);
OMEGACUTpre = [0,2];
OMEGACUT = OMEGACUTpre*OmegaFactor;
%
bandplot(F0CAR,OMEGACUTpre*OmegaFactor,'ylabel','Frequency (Hz)','ax',FigTest.axes(2));
