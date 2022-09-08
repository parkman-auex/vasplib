%% GreenFunction in HR
% BHZ model
%% BHZ model
% Useful matrix

%% useful tool
s_0   = pauli_matric(0)  ;  s_x = pauli_matric(1);  s_y =  pauli_matric(2) ;  s_z = pauli_matric(3);
sigma_0 = pauli_matric(0);sigma_x =  pauli_matric(1);sigma_y =  pauli_matric(2);sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0);  tau_x =  pauli_matric(1);  tau_y =  pauli_matric(2);  tau_z = pauli_matric(3);
% k·p model

syms C0 C1 C2 real;
syms M0 M1 M2 real;
syms A B real;
syms k_x k_y k_z real;
%
M       = M0-M2*(k_x^2+k_y^2)-M1*(k_z^2);
E0k     = C0+C2*(k_x^2+k_y^2)+C1*(k_z^2);
k_plus  = k_x + 1i* k_y;
k_minus = k_x - 1i* k_y;
%
BHZ = HK(4,2)
%
BHZ = BHZ ...
    +Term(A*k_x ,sigma_z*tau_x )...
    +Term(A*k_y ,sigma_0*tau_y )...
    +Term(E0k   ,sigma_0*tau_0 )...
    +Term(M     ,sigma_0*tau_z )...
    +Term(B*(k_x^2-k_y^2),sigma_x*tau_x )...
    +Term(2*B*(k_x*k_y),sigma_y*tau_x )...
    ;
disp(BHZ);
% Bulk band
% para

% unit eV 
C0     = -0.0   ;
M0     = -1   ;
%M0     = -0.25   ;
% unit eV Ang
A      =  1      ;
%B      =  4.1      ;
% unit eV Ang^2
C1     =  0     ;
C2     =  0     ;
M1     =  -0.5      ;
M2     =  -0.5  ;
B     =  0    ;
% unit    Ang
a      =  1        ;
b      =  1       ;
c      =  1        ;
% bandplot

BHZ_n = BHZ.Subsall();
BHZ_n = BHZ_n <'POSCAR';
BHZ_n = BHZ_n <'KPOINTS';
EIGENCAR = BHZ_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-6,6]);
% kp2TB

BHZ = BHZ <'POSCAR';
BHZ_TB= BHZ.kp2TB();
C1     =  0     ;
BHZ_TB_n = BHZ_TB.Subsall();
BHZ_TB_n = BHZ_TB_n <'POSCAR';
BHZ_TB_n = BHZ_TB_n <'KPOINTS';
BHZ_TB_n.bandplot([-3,3])
%% Slab-band

repeatnum   = 20;
fin_dir     =  2;
[EIGENCAR_slab,klist_l,kpoints_l,kpoints_name] = BHZ_TB_n.slab(repeatnum,fin_dir,'KPOINTS_1D_X');
% plot
bandplot(EIGENCAR_slab,[-3,3],klist_l,kpoints_l,kpoints_name,'title','BHZ-slab','Color','g');
%% Surf Green
% Three ways to gen surf spectrum GREENCAR
% Tmatrix iter

w_range = [-2,2,500];
fin_dir = 2 ;
KPOINTS_surf = 'KPOINTS_1D_X';
principle_layer = 1;
eta =0.01;
fermi = 0;
mode = 'Tmatrix_iter';
tic;
[DOSCAR_l,~,~,w_list,klist_l,kpoints_l,kpoints_name] = ...
   BHZ_TB_n.surf(w_range,fin_dir,KPOINTS_surf,principle_layer,eta,fermi,mode);
toc;
heatplot(DOSCAR_l,w_list,klist_l,kpoints_l,kpoints_name);
% Tmatrix inv

w_range = [-2,2,500];
fin_dir = 2 ;
KPOINTS_surf = 'KPOINTS_1D_X';
principle_layer = 1;
eta =0.01;
fermi = 0;
mode = 'Tmatrix_inv';
 warning off;

tic;
[DOSCAR_l,~,~,w_list,klist_l,kpoints_l,kpoints_name] = ...
   BHZ_TB_n.surf(w_range,fin_dir,KPOINTS_surf,principle_layer,eta,fermi,mode);
toc;
heatplot(DOSCAR_l,w_list,klist_l,kpoints_l,kpoints_name);
%% 
% 这种方法对模型要求高, 并非总是有效，有时候会比第三种方法快
% Green iter(default)

% mode = 'GW_iter';
% warning on;
w_range = [-2,2,500];
eta =0.01;
tic;
[DOSCAR_l,~,~,w_list,klist_l,kpoints_l,kpoints_name] = ...
   BHZ_TB_n.surf(w_range,fin_dir,KPOINTS_surf,principle_layer,eta);
toc;
heatplot(DOSCAR_l,w_list,klist_l,kpoints_l,kpoints_name);
%% 
% 这种默认方法最快
% Three ways to gen surf fermi arc GREENCAR
% break chiral

C1 = 0.1;
BHZ_TB_n = BHZ_TB.Subsall();
EIGENCAR = BHZ_TB_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-6,6]);

ax = gca;
chart = ax.Children(1);
datatip(chart,1.597,0.2213);
% Tmatrix iter

w_arc = 0.;
fin_dir = 2 ;
kmesh = [50 50];
kfermiarc = [-0.1, 0,-0.5;... %begin kpoints
              0.2 , 0,   0;... %  first vector
              0   , 0,   1;... %  second vector
             ];
principle_layer = 1;
eta =0.01;
fermi = 0.231945;
mode = 'Tmatrix_iter';
%
tic;
[DOSCAR_l,~,~,klist1,klist2] = ...
   BHZ_TB_n.fermiarc(w_arc,fin_dir,kmesh,kfermiarc,principle_layer,eta,fermi,mode);
toc;
heatplot(DOSCAR_l,klist2,klist1,'k_x','k_z',turbo,'arc');
% Tmatrix inv

fermi = 0.231945;
mode = 'Tmatrix_inv';
%
warning off;
tic;
[DOSCAR_l,~,~,klist1,klist2] = ...
   BHZ_TB_n.fermiarc(w_arc,fin_dir,kmesh,kfermiarc,principle_layer,eta,fermi,mode);
toc;
heatplot(DOSCAR_l,klist2,klist1,'k_x','k_z',turbo,'arc');
% GreenFunction iter

fermi = 0.231945;
%
warning off;
tic;
[DOSCAR_l,~,~,klist1,klist2] = ...
   BHZ_TB_n.fermiarc(w_arc,fin_dir,kmesh,kfermiarc,principle_layer,eta,fermi);
toc;
heatplot(DOSCAR_l,klist2,klist1,'k_x','k_z',turbo,'arc');
% Three ways to gen surf fermi arc*(3D) GREENCAR

C1 = 0.1;
BHZ_TB_n = BHZ_TB.Subsall();
w_range = [-0.5,0.5,50];
fin_dir = 2 ;
kmesh = [50 60];
kfermiarc = [-0.1, 0,-0.5;... %begin kpoints
              0.2 , 0,   0;... %  first vector
              0   , 0,   1;... %  second vector
             ];
principle_layer = 1;
eta =0.01;
fermi = 0.231945;
[DOSCAR_l,~,~,klist1,klist2] = ...
   BHZ_TB_n.fermiarc3D(w_range,fin_dir,kmesh,kfermiarc,principle_layer=1,eta=0.01,fermi=0.231945);
%%