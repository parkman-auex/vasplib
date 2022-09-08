%% HRclass Tutorial
% 基于HRclass实现 WSM model      
% 
% 2021.08
%% 
% * Author: parkman
% * Email：parkman@buaa.edu.cn
%% Model
% Prepare

% useful tool
s_0   = pauli_matric(0)  ;  s_x = pauli_matric(1);  s_y =  pauli_matric(2) ;  s_z = pauli_matric(3);
sigma_0 = pauli_matric(0);sigma_x =  pauli_matric(1);sigma_y =  pauli_matric(2);sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0);  tau_x =  pauli_matric(1);  tau_y =  pauli_matric(2);  tau_z = pauli_matric(3);

syms C0 C1 C2 real;
syms M0 M1 M2 real;
syms A0 real;
syms k_x k_y k_z real;
% kp

M       = M0-M2*(k_x^2+k_y^2)-M1*(k_z^2);
E0k     = C0+C2*(k_x^2+k_y^2)+C1*(k_z^2);
A       = A0;%+A2*(k_x^2+k_y^2)+A1*(k_z^2)
k_plus  = k_x + 1i* k_y;
k_minus = k_x - 1i* k_y;
%
QWZ_3D = HK(2,2);
QWZ_3D = QWZ_3D ...
    +Term(A*k_x ,tau_x )...
    +Term(A*k_y ,tau_y )...
    +Term(E0k   ,tau_0 )...
    +Term(M     ,tau_z );...
QWZ_3D = QWZ_3D.Subsall('sym');
QWZ_3D = QWZ_3D <'POSCAR_4';
%%
QWZ_3D_kp = QWZ_3D.sym()
%QWZ_3D_kp_pauli = QWZ_3D.pauliDecomposition
% kp2TB

QWZ_3D_TB= QWZ_3D.kp2TB();
QWZ_3D_TB.list();
% Para

% unit eV 
M0     = 1   ;
%M0     = -0.25   ;
% unit eV Ang
A0      =  1      ;
%B      =  4.1      ;
% unit eV Ang^2
C0     =  0;
C1     =  0.1     ;
C2     =  0     ;
M1     =  0.5      ;
M2     =  0.5  ;
% unit    Ang
a      =  1        ;
b      =  1       ;
c      =  1        ;
% NumTB

QWZ_3D_TB_n = QWZ_3D_TB.Subsall();
QWZ_3D_TB_n = QWZ_3D_TB_n <'KPOINTS_4';
%% Bulk
% bandstructrue

EIGENCAR = QWZ_3D_TB_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-2,2],'title',"WSM-TETRA",'POSCAR','POSCAR_4','KPOINTS','KPOINTS_4');
% 3d Band

[EIGENCAR_3D,klist1,klist2] = QWZ_3D_TB_n.EIGENCAR_gen_3D([101,101],[-0.5,-0.0,-0.5;1,0,0;0 0 1],'fin_dir',2);
[fig,ax] = QWZ_3D_TB_n.BZplot(QWZ_3D_TB_n.Rm,'mode','3D');
vasplib.bandplot_3d(EIGENCAR_3D,klist1,klist2,'ax',ax,'xlabel','k_x','ylabel','k_z','cmap',hsv);view(70,9)
%% Surf
% slabband

QWZ_3D_TB_n_slab = QWZ_3D_TB_n.cut_piece(50,2);
QWZ_3D_TB_n_slab = QWZ_3D_TB_n_slab<'KPOINTS_slab';
ProjectionStruct.discrimination = 0.1;
ProjectionStruct.center = [0.5,0.5,0.5];
ProjectionStruct.orientation = 2;
ProjectionStruct.sign = true;
%%
[EIGENCAR_slab,~,WEIGHTCAR_slab]= QWZ_3D_TB_n_slab.EIGENCAR_gen('LWAVE',false,'WEIGHTCAR','true','ProjectionMethod','slab','ProjectionStruct',ProjectionStruct);
%%
%gnuplot_pm3d = ColorMap.gnuplot.pm3d;
%ColorMap.Matplotlib();
coolwarm = ColorMap.Matplotlib('coolwarm');
%red_gray_blue = ColorMap.customcolormap([0  .5  1], {'#68011d','#808080','#062e61'});
vasplib.pbandplot(WEIGHTCAR_slab,EIGENCAR_slab,'POSCAR','POSCAR_4','KPOINTS','KPOINTS_slab','cmap',coolwarm );
% slabband_3d

ProjectionStruct.sign = true;
[EIGENCAR_3D,klist1,klist2,WEIGHTCAR_3D] = QWZ_3D_TB_n_slab.EIGENCAR_gen_3D([101,101],[-0.1,-0.0,-0.5;0.2,0,0;0 0 1],'fin_dir',2,...
'WEIGHTCAR',true,'ProjectionMethod','slab','ProjectionStruct',ProjectionStruct);
%%
%[fig,ax] = QWZ_3D_TB_n.BZplot(QWZ_3D_TB_n.Rm,'mode','2D');
gnuplot_pm3d = ColorMap.gnuplot.pm3d;
vasplib.bandplot_3d(EIGENCAR_3D,klist1,klist2,'WEIGHTCAR',WEIGHTCAR_3D,'xlabel','k_x','ylabel','k_z','cmap',coolwarm);view(70,9)
% slice

Earc = 0.2;
infinite_small = 1e-2;
EIGENCAR_3D_slice = EIGENCAR_3D;
WEIGHTCAR_3D_slice = WEIGHTCAR_3D;
ChooseXY=find(EIGENCAR_3D_slice >Earc +infinite_small | EIGENCAR_3D <Earc -infinite_small  );
EIGENCAR_3D_slice(ChooseXY) =0;
WEIGHTCAR_3D_slice(ChooseXY) =nan;
vasplib.bandplot_3d(EIGENCAR_3D_slice,klist1,klist2,'WEIGHTCAR',WEIGHTCAR_3D_slice,'xlabel','k_x','ylabel','k_z','title','E = 0.2 eV slice','cmap',coolwarm);