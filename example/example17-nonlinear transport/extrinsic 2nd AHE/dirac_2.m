%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%%
syms t1 t2 v1 v2 eta1 eta2 m E1 E2;
syms gamma;
syms k_x k_y real;
H_bare = HK(2,1);
%%
H_d1 = H_bare ...
    + Term( t1*k_x + E1,    sigma_0)...
    + Term( v1*k_y,         sigma_x)...
    + Term( v1*eta1*k_x,    sigma_y)...
    + Term( m/2,            sigma_z);
H_d2 = H_bare ...
    + Term( t2*k_x + E2,    sigma_0)...
    + Term( v2*k_y,         sigma_x)...
    + Term( v2*eta2*k_x,    sigma_y)...
    + Term( m/2,            sigma_z);

anti_I4 = kron([0 1;1 0],[0 1;1 0]);
H_dirac = [H_d1,H_d2] + Term(gamma,anti_I4);
%%
t1 = 0.5;
t2 = 0;
v1 = 1;
v2 = 1;
eta1 = -1;
eta2 = 1;
m = 0.4;
E1 = 0.05;
E2 =-0.05;
gamma = 0.01;
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