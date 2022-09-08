%% Quadratic and cubic nodal lines stabilized by crystalline symmetry
% https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.121106 (5)
% Heff(k) = a * k_-^2+ H.c.,
%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
sigma_plus = double(sigma_x) + 1i*double(sigma_y);
sigma_minus = double(sigma_x) - 1i*double(sigma_y);
%% 
syms a u1 u2 real;
syms k_x k_y k_z d real ;
k_minus = (k_x - 1i*k_y) ;
k_plus = (k_x + 1i*k_y) ;
H_174 = HK(2,2);
%%
H_174 = H_174 ...
    + Term( a*k_minus^2 ,    sigma_plus)...
    + Term( a*k_plus^2  ,    sigma_minus)...
    ;
% H_174  =H_174  ...
%     + Trig( u1*(cos(d*k_z) )         ,    sigma_x)...
%     + Trig( u2*(sin(d*k_z) )         ,    sigma_y)...
%     ;
%%
H_174 = H_174  < 'POSCAR';
H_174 = H_174  < 'KPOINTS';
d =3.3;
a = 1;
% u1 =0.233;
% u2 =0.002;
u1 =0;
u2= 0;
tolerance = 0.01;
H_174n = H_174.Subsall();
%
% EIGENCAR = H_174n.EIGENCAR_gen();
% bandplot(EIGENCAR,[-1,1]);
H_174_TB = H_174.kp2TB();
% H_174_TB = H_174_TB <'POSCAR';
H_174_TB = H_174_TB <'KPOINTS';
H_174_TBn = H_174_TB.Subsall();
EIGENCAR_TB = H_174_TBn.EIGENCAR_gen();
bandplot(EIGENCAR_TB,[-1,1]);
%%
% kzf = 0:0.01:1;
% H_174n  = H_174n < 'KPOINTS_findnode';
% [k_list_cart,klist_f,gaplist]=H_174n.findnodes2(kzf,1,tolerance);
%%
GammaOper = double(sigma_z);
kloop = H_174_TBn.kloop_gen([0 0 0;0 0 1;0.1 0 101],'kplane');
kpoints_r = [0 0 0.3];
H_174n.TopoCharge(GammaOper,kloop,kpoints_r);

