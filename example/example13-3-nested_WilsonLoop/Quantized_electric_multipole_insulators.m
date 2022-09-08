clear;					
%% Quantized electric multipole insulators 
s_0   = pauli_matric(0)  ;  s_x = pauli_matric(1);  s_y =  pauli_matric(2) ;  s_z = pauli_matric(3);
sigma_0 = pauli_matric(0);sigma_x =  pauli_matric(1);sigma_y =  pauli_matric(2);sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0);  tau_x =  pauli_matric(1);  tau_y =  pauli_matric(2);  tau_z = pauli_matric(3);
Gamma_0 = tau_z*sigma_0;
Gamma_1 = -1*tau_y*sigma_x;
Gamma_2 = -1*tau_y*sigma_y;
Gamma_3 = -1*tau_y*sigma_z;
Gamma_4 = -1*tau_x*sigma_0;
syms gamma lambda delta real;
syms k_x k_y k_z real ;
%
%% 
% 
QTI = Htrig(4);
%
QTI = QTI ...
    +Trig(gamma + lambda*cos(k_x) ,Gamma_4 )...
    +Trig(lambda*sin(k_x) ,Gamma_3 )...
    +Trig(gamma + lambda*cos(k_y) ,Gamma_2 )...
    +Trig(lambda*sin(k_y) ,Gamma_1 )...s_z
    +Trig(delta,Gamma_0)...
    ;
QTI_sym = sym(QTI)
QTI = QTI.Subsall('sym');
QTI = QTI <'POSCAR_4';
QTI = QTI <'KPOINTS_4';
QTI_hr = QTI.Htrig2HR();
%%
delta = 0;
lambda = 1;
gamma = 0.5;
% bandplot
QTI_n = QTI.Subsall();
QTI_hr_n = QTI_hr.Subsall();
% QTI_n = QTI_n <'POSCAR_4';
% QTI_n = QTI_n <'KPOINTS_4';
%%
Fig = Figs(1,2);
bandplot(QTI_n.EIGENCAR_gen(),'title',"QTI-TETRA-Htrig",'KPOINTS','KPOINTS_4','ax',Fig.axes(1));
bandplot(QTI_hr_n.EIGENCAR_gen(),[-3,3],'title',"QTI-TETRA-Hr",'KPOINTS','KPOINTS_4','ax',Fig.axes(2));
%%
[BFCAR,~,klist_l] = QTI_n.WilsonLoop('knum_evol',1);
vasplib.WilsonLoopPlot(BFCAR,klist_l);
%%
[BFCAR,~,klist_l] = QTI_hr_n.WilsonLoop('knum_evol',1);
vasplib.WilsonLoopPlot(BFCAR,klist_l);
%%
[nested_BFCAR,~,nested_klist_l] = QTI_n.nested_WilsonLoop('knum_evol',2,'knum_int',51);
 vasplib.WilsonLoopPlot(nested_BFCAR,nested_klist_l);
%%
[nested_BFCAR,~,nested_klist_l] = QTI_hr_n.nested_WilsonLoop('knum_evol',2,'knum_int',101);
vasplib.WilsonLoopPlot(nested_BFCAR,nested_klist_l);
%%
QTI_n2 = QTI_n;
QTI_n2.orbL = zeros(4,3);
[nested_BFCAR,~,nested_klist_l] = QTI_n2.nested_WilsonLoop();
vasplib.WilsonLoopPlot(nested_BFCAR,nested_klist_l);