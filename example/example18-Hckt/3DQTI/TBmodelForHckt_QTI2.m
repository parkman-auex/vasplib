%% armchair 2
Type = 'list';%listmat
syms v w u real;
QTI = HR.from_POSCAR_SE('POSCAR','level_cut',5,'chiral',true,'onsite',false,'per_dir',[1,1,0],'WAN_NUM',4,'r_max',2,'Type',Type);
QTI = subs(QTI,[v u w]);
QTI = QTI<'POSCAR';

selectL = ((QTI.vectorL(:,4) == 3)&(QTI.vectorL(:,5) == 1) ) ...
    | ((QTI.vectorL(:,4) == 1)&(QTI.vectorL(:,5) == 3) );
QTI.HcoeL(selectL) = -QTI.HcoeL(selectL);

%QTI = QTI.autohermi();

HcktTB = QTI;
%%

%% Para

v = 1; % transverse bond
w = 1.5;
C_0 = 1;
u = 1;
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
