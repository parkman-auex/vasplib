%% BHZ kp band and berry phase
BHZ_kp;
% unit eV 
C0     = -0.0   ;
M0     = -1   ;
%M0     = -0.25   ;
% unit eV Ang
A0      =  1      ;
%B      =  4.1      ;
% unit eV Ang^2
C1     =  0     ;
C2     =  0     ;
M1     =  -0.5      ;
M2     =  -0.5  ;
B     =  0.5    ;
% unit    Ang
a      =  1        ;
b      =  1       ;
c      =  1        ;
% bandplot

BHZ_n = BHZ.Subsall();
BHZ_n = BHZ_n <'POSCAR_BHZ';
BHZ_n = BHZ_n <'KPOINTS_BHZ';
EIGENCAR = BHZ_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-6,6], ...
    'title','BHZ-HODSM', ...
    'POSCAR','POSCAR_BHZ', ...
    'KPOINTS','KPOINTS_BHZ');
%%
[BFCAR,~,klist_l] = BHZ_n.WilsonLoop('kevolution',[2,0,0]);
[fig,ax] = vasplib.WilsonLoopPlot(BFCAR,klist_l);