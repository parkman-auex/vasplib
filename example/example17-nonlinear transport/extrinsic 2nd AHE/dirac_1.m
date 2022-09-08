%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%%
syms t v eta m alpha;
syms k_x k_y k_z real;
H_dirac = HK(2,2);
%%
H_dirac = H_dirac ...
    + Term( t*k_x,                   sigma_0)...
    + Term( v*k_y,                   sigma_x)...
    + Term( v*eta*k_x,               sigma_y)...
    + Term( m/2-alpha*(k_x^2+k_y^2), sigma_z)...
    ;

H_dirac = H_dirac < 'KPOINTS';
H_dirac = H_dirac.input_Rm();
%%
t = 0.5;
v = 1;
eta = -1;
m = 0.2;
alpha = 1;
H_dirac_n = H_dirac.Subsall();
%%
% EIGENCAR = H_dirac_n.EIGENCAR_gen();
% bandplot(EIGENCAR);
%%
nk = [80,80,1];
vk = [0.2 0   0;
      0   0.2 0;
      0   0   1];
[~,~,~,~,kstruct] = kmesh_gen(H_dirac_n,"nk",nk,'vk',vk,...
    "mode","center","edge","full");
%%
D_bd = transport.SOAHC_ext(H_dirac_n,kstruct,"xz",...
    'Ef_range',[-1,1],'T',50);