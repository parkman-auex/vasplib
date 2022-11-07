sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
%% 
syms t;
t_so = 0.1*t;
syms k_x k_y k_z real ;
a = 1.42;
f_x = t*(1+2*cos(3*a*k_y/2)*cos(sqrt(3)*a*k_x/2));
f_y = t*2*sin(3*a*k_y/2)*cos(sqrt(3)*a*k_x/2);
f_so= -2*t_so*(sin(sqrt(3)*a*k_x)-2*cos(3*a*k_y/2)*sin(sqrt(3)*a*k_x/2));
%%
Hal = Htrig(2);
Hal = Hal ...
    + Trig( f_x,    sigma_x)...
    + Trig( f_y,    sigma_y)...
    + Trig( f_so,   sigma_z);
Hal = Hal < 'POSCAR';
Hal = Hal < 'KPOINTS';
%%
t=0.8;
Hal_n = Hal.Subsall();
EIGENCAR = Hal_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-3,3],'title',"Haldane");
%% vasplin inner function
chern0 = Hal_n.Chern()
chern1 = Hal_n.Chern('BAND_index',1)
chern2 = Hal_n.Chern('BAND_index',2)
[BFCAR,~,klist_l] = Hal_n.WilsonLoop('knum_evol',100,'cartesian',true);
vasplib.WilsonLoopPlot(BFCAR,klist_l)