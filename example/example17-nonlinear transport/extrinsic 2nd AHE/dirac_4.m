%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%%
syms t1 t2 v1 v2 eta1 eta2 m E1 E2;
syms gamma vm vn;
syms k_x k_y real;
k_x_1 = k_x + 0.1*pi;
k_x_2 = k_x + 0.15*pi;
H0 = HK(4,1);
%%
H_d1 = H0 ...
    + Term( t1*k_x_1 + E1,  sigma_0*sigma_0)...
    + Term( v1*k_y,         sigma_0*sigma_x)...
    + Term( v1*eta1*k_x_1,  sigma_0*sigma_y)...
    + Term( m/2,            sigma_0*sigma_z);
%%
H_d2 = H0 ...
    + Term( t2*k_x_2 + E2,  sigma_0*sigma_0)...
    + Term( v2*k_y,         sigma_0*sigma_x)...
    + Term( v2*eta2*k_x_2,  sigma_0*sigma_y)...
    + Term( m/2,            sigma_0*sigma_z);
%%
H_p = H0 ...
    + Term( vm*k_x,         sigma_x*sigma_z)...
    + Term(-1i*vn*k_y,      sigma_x*sigma_0);
%%
H_dirac = [H_d1+H_p, H_d2+H_p] + Term(gamma, sigma_x*sigma_x*sigma_x);
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
vm = 0.05;
vn = 0;
H_dirac_n = H_dirac.Subsall();
%%
% EIGENCAR = H_dirac_n.EIGENCAR_gen();
% bandplot(EIGENCAR);
%%
nk = [100,100,1];
[~,~,~,~,kstruct] = kmesh_gen(H_dirac_n,[0 0 0; 0.1 0 0; 0 0.1 0; 0 0 1],...
    "nk",nk,"mode","center","edge","full");
%%
D_bd = transport.SOAHC_ext(H_dirac_n,kstruct,"xz",...
    'Ef_range',[-1,1],'T',50);