%% 
clear;
% 10.1103/PhysRevLett.121.266601
%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%% dour dirac hk 
syms t1 t2 v1 v2 eta1 eta2 m E1 E2 real;
syms gamma v_x v_y alpha real;
syms k_x k_y real;
k_x_1 = k_x + 0.1*pi;
k_x_2 = k_x + 0.15*pi;
k_2 = (k_x +k_y)^2;
H0 = HK(4,1);
%%
H_d1 = H0 ...
    + Term( t1*k_x_1 + E1,  sigma_0*sigma_0)...
    + Term( v1*k_y,         sigma_0*sigma_x)...
    + Term( v1*eta1*k_x_1,  sigma_0*sigma_y)...
    + Term( m/2,  sigma_0*sigma_z);
%%
H_d2 = H0 ...
    + Term( t2*k_x_2 + E2,  sigma_0*sigma_0)...
    + Term( v2*k_y,         sigma_0*sigma_x)...
    + Term( v2*eta2*k_x_2,  sigma_0*sigma_y)...
    + Term( m/2,  sigma_0*sigma_z);
%%
H_p = H0 ...
    + Term( v_x*k_x,         sigma_x*sigma_z)...
    + Term( -v_y*k_y,        sigma_y*sigma_0);
%%
H_all = [H_d1+H_p, H_d2+H_p] + Term(gamma, sigma_x*sigma_x*sigma_x);
H_all_sym = sym(H_all); % better in realtime script
H_all = H_all <'POSCAR_8band';
%%
t1 = 1.5;
t2 = 1.5;
v1 = 2;
v2 = 2;
eta1 = -1;
eta2 = 1;
m = 0.1;
E1 = 0.02;
E2 =-0.08;
gamma = 0.05;
v_x = 0.0;
v_y = 0;
alpha = 0 ;
H_all_n = H_all.Subsall();
H_all_n = H_all_n <'KPOINTS_8band';
EIGENCAR = H_all_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-0.5,0.5],'POSCAR','POSCAR_8band','KPOINTS','KPOINTS_8band');
%%
sizes = 101;
[BCCAR,Grid,~,klist] = H_all_n.BC_2D( ...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-1,-1,0],'kdir1',[2,0,0],'kdir2',[0,2,0],...
    'plot',false,'oneshot',true);
[fig,ax] = vasplib.BCplot2D(BCCAR,Grid,H_all_n.Rm,'BZ',true);
fprintf('Chern Number is %7.5f\n',sum(BCCAR,'all')/(2*pi));
%%
% figure(2)
% Ef_range = linspace(-0.5,0.5,101);
% [D_xz,omega_xy] = BCD_kubo(H_all_n,Ef_range);
% plot(D_xz,Ef_range)