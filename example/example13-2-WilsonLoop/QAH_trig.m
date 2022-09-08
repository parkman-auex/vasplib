%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%% 
syms t1 t2 t3 m;
syms k_x k_y k_z real ;
d1 = 2*t1*cos(k_x);
d2 = 2*t1*cos(k_y);
d3 = m + 2*(t3*sin(k_x) + t3*sin(k_y) + t2*cos(k_x+k_y));
%%
QAH = Htrig(2);
QAH = QAH ...
    + Trig(d1,  sigma_x)...
    + Trig(d2,  sigma_y)...
    + Trig(d3,  sigma_z);
%%
t1 = 1.0;
t2 = 1.0;
t3 = 0.5;
m = -1.0;
QAH_n = QAH.Subsall();
chern0 = chern_number(QAH_n);
chern1 = chern_number(QAH_n,1);
chern2 = chern_number(QAH_n,2);
wilson_loop(QAH_n,"kx");
%% vasplin inner function
chern0 = QAH_n.Chern()
chern1 = QAH_n.Chern('BAND_index',1)
chern2 = QAH_n.Chern('BAND_index',2)
[BFCAR,~,klist_l] = QAH_n.WilsonLoop();
vasplib.WilsonLoopPlot(BFCAR,klist_l)
