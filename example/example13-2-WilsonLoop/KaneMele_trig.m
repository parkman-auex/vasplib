%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%% 
syms t;
t_so = 0.1*t;
lambda = 0.2*t;
syms k_x k_y k_z real ;
a = 1.42;
f_x = t*(1+2*cos(3*a*k_y/2)*cos(sqrt(3)*a*k_x/2));
f_y = t*2*sin(3*a*k_y/2)*cos(sqrt(3)*a*k_x/2);
f_so= -2*t_so*(sin(sqrt(3)*a*k_x)-2*cos(3*a*k_y/2)*sin(sqrt(3)*a*k_x/2));
%%
KM = Htrig(4);
KM = KM ...
    + Trig( f_x,    sigma_x*tau_0)...
    + Trig( f_y,    sigma_y*tau_0)...
    + Trig( f_so,   sigma_z*tau_z)...
    + Trig( lambda, sigma_0*tau_y);
KM = KM < 'POSCAR';
KM = KM < 'KPOINTS';
%%
t=0.8;
KM_n = KM.Subsall();
% EIGENCAR = KM_n.EIGENCAR_gen();
% bandplot(EIGENCAR,[-3,3],"KaneMele");
chern1 = chern_number(KM_n);
chern2 = chern_number(KM_n,1:2);
wilson_loop(KM_n,"kx");
%% vasplin inner function
chern0 = KM_n.Chern()
chern1 = KM_n.Chern('BAND_index',1)
chern2 = KM_n.Chern('BAND_index',2)
% bugs here
[BFCAR,~,klist_l] = KM_n.WilsonLoop('cartesian',true);
vasplib.WilsonLoopPlot(BFCAR,klist_l)