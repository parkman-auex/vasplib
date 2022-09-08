%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%% 
syms tt t_sp ee g lambda B_x B_y;
syms k_x k_y k_z real ;
f_2 = ee - tt*(cos(k_x)+cos(k_y));
%%
BHZ = Htrig(4);
BHZ = BHZ ...
    + Trig( f_2,    sigma_z*tau_0)...
    + Trig( pi/2 *2*t_sp*sin(k_x), sigma_x*tau_z)...
    + Trig( pi/2 *2*t_sp*sin(k_y), sigma_y*tau_0)...
    + Trig( pi/2 *g*lambda*B_x*sin(k_z)*(cos(k_x)-cos(k_y)), sigma_x*tau_x)...
    + Trig( pi/2 *g*lambda*B_y*sin(k_z)*sin(k_x)*sin(k_y), sigma_y*tau_x)...
    ;
BHZ = BHZ < 'POSCAR';
BHZ = BHZ < 'KPOINTS';
%%
tt = 1; 
ee = 1;
t_sp = 0.3;
g = 1;
lambda = 1; 
B_x =0.15;
B_y =0.15;
BHZ_n = BHZ.Subsall();
% EIGENCAR = BHZ_n.EIGENCAR_gen();
% bandplot(EIGENCAR,[-3,3],"BHZ");
% !!! not support the spin Chern number of BHZ model
%wilson_loop(BHZ_n,"kx");
[BFCAR,~,klist_l] = BHZ_n.WilsonLoop('kstart',[0,0,0.1]);
[fig,ax] = vasplib.WilsonLoopPlot(BFCAR,klist_l);