%% 创建一个SSHTB模型

SSH = HR.from_POSCAR_SE('POSCAR',...
'Type','list',...
'onsite',false,...
'Chiral',true,...
'search_range' ,[1 0 0],...
'r_max',5,...
'level_cut',2);
SSH.list;
%%
syms C_v C_w Minus real;
SSH = SSH.subs([C_v C_w]);
Gauge = 1;
switch Gauge
    case 1
        SSH.HcoeL(SSH.HcoeL ==C_v) = Minus*C_v;

end
HcktTB = SSH;
%
%% Para
Minus = -1;
C_inter = 1;
C_v = 1;
C_w = 2;
C_0 = 3;
%% prepare
Hcktpre = HcktTB.HRforHckt();
%% Add Unit
L_0 = 1e-6;
MU = 100;
try
    switch magnitude
        case 'p'
            Cfactor = MU*1E-12;% 100 pF
        case 'n'
            Cfactor = MU*1E-09;% 100 nF
        case 'u'
            Cfactor = MU*1E-06;% 100 uF
        case 'm'
            Cfactor = MU*1E-03;% 100 mF
    end
catch
    Cfactor = 100*1E-12;% 100 pF
end
OmegaFactor = (Cfactor*L_0*100)^(-1/2);

%%
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