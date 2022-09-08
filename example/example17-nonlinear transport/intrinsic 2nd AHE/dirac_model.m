%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%%
syms w vm vn delta;
syms k_x k_y k_z real;
H_dirac = HK(4,1);
%%
H_dirac = H_dirac ...
    + Term( w * k_x,       tau_0 * sigma_0)...
    + Term( vm * k_x,      tau_x * sigma_0)...
    + Term( vn * k_y,      tau_y * sigma_y)...
    + Term( delta,         tau_z * sigma_0);
%%
w = 0.4;
vm = 1;
vn = 1;
delta = 40 * 1e-3;
H_dirac_n = H_dirac.Subsall();
%%
[~,~,~,~,kstruct] = kmesh_gen(H_dirac_n,...
    [0 0 0; 0.05 0 0; 0 0.05 0; 0 0 1],...
    "nk",[120 200 1], "mode","center", "edge","full");
%%
Chi_abc_HK = transport.SOAHC_int(H_dirac_n,kstruct,"xyy",...
    'Ef_range',[-0.3,0.3],...
    'plotBZ_bands',1:2);