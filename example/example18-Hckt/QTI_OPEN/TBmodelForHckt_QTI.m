%% armchair 2

sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

QTI = HR(4);
QTI = QTI<'POSCAR_4orb';
QTI = QTI<'KPOINTS';

% orb_order
% 1 t1 2 t1
% m1  t2
% 3 t1 4 t1
% m2  t3

V_N21 = [0 0 0;1 0 0];
V_N12 = [0 0 0;-1 0 0];
V_N24 = [0 0 0; 0 1 0];
V_N42 = [0 0 0; 0 -1 0];

syms t1 t2 t3 m1 m2 tm1 tm2 l1 l2 real;

QTI = QTI.set_hop(...
     t1,1,2,[0 0 0],'symadd');
QTI = QTI.set_hop(...
     t1,2,1,[0 0 0],'symadd');
QTI = QTI.set_hop(...
     t2,1,2,[-1 0 0],'symadd');
QTI = QTI.set_hop(...
     t2,2,1,[1 0 0],'symadd');

QTI = QTI.set_hop(...
     tm1,3,4,[0 0 0],'symadd');
QTI = QTI.set_hop(...
     tm1,4,3,[0 0 0],'symadd');
QTI = QTI.set_hop(...
     tm2,3,4,[-1 0 0],'symadd');
QTI = QTI.set_hop(...
     tm2,4,3,[1 0 0],'symadd');


QTI = QTI.set_hop(...
     l1,2,4,[0 0 0],'symadd');
QTI = QTI.set_hop(...
     l1,4,2,[0 0 0],'symadd');
QTI = QTI.set_hop(...
     l2,2,4,[0 1 0],'symadd');
QTI = QTI.set_hop(...
     l2,4,2,[0 -1 0],'symadd');


QTI = QTI.set_hop(...
     m1,1,3,[0 0 0],'symadd');
QTI = QTI.set_hop(...
     m1,3,1,[0 0 0],'symadd');
QTI = QTI.set_hop(...
     m2,1,3,[0 1 0],'symadd');
QTI = QTI.set_hop(...
     m2,3,1,[0 -1 0],'symadd');


%QTI = QTI.autohermi();

HcktTB = QTI;

 %% Para

t1 = 1; % transverse bond
t2 = 1.5;
tm1 = 1; % transverse bond
tm2 = 1.5;


m1 = -1;
m2 = -1.5;
l1 = 1;
l2 = 1.5;

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
    Cfactor = 100*1E-6;% 100 pF
end
OmegaFactor = (Cfactor*L_0*100)^(-1/2);

C_0 = 1;
% drawnow;
FigTest = Figs(1,2);
bandplot(HcktTB.Subsall().EIGENCAR_gen('convention','I'),...
    [-5,5],'Color','b','ax',FigTest.axes(1));
%
Hcktpre_n = Hcktpre.Subsall();
%
EIGENCAR = Hcktpre_n.EIGENCAR_gen();
F0CAR = (EIGENCAR*L_0*Cfactor).^(-1/2)./(2*pi);
OMEGACUTpre = [0.5,1];
OMEGACUT = OMEGACUTpre*OmegaFactor;
%
bandplot(F0CAR,OMEGACUTpre*OmegaFactor,'ylabel','Frequency (Hz)','ax',FigTest.axes(2));
