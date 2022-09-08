% BBH
% ladimir A. Benalcazar, B. Andrei Bernevig, and Taylor L. Hughes, "Quantized electric multipole insulators", Science Volume 357:61-66 (2017).
%% Quantized electric multipole insulators 
s_0   = pauli_matric(0)  ;  s_x = pauli_matric(1);  s_y =  pauli_matric(2) ;  s_z = pauli_matric(3);
sigma_0 = pauli_matric(0);sigma_x =  pauli_matric(1);sigma_y =  pauli_matric(2);sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0);  tau_x =  pauli_matric(1);  tau_y =  pauli_matric(2);  tau_z = pauli_matric(3);
Gamma_0 = tau_z*sigma_0;
Gamma_1 = -1*tau_y*sigma_x;
Gamma_2 = -1*tau_y*sigma_y;
Gamma_3 = -1*tau_y*sigma_z;
Gamma_4 = -1*tau_x*sigma_0;
syms gamma lambda Delta_x Delta_y delta_x delta_y real;
syms k_x k_y k_z real ;
%
gamma_x = gamma;
gamma_y = gamma;
lambda_x = lambda;
lambda_y = lambda;
%%
H = zeros(4,'sym');
H(1, 1) = Delta_x;
H(2, 2) = Delta_y;
H(3, 3) = Delta_y;
H(4, 4) = Delta_x;
H(1, 3) = gamma_x+lambda_x*exp(1i*k_x);
H(2, 4) = gamma_x+lambda_x*exp(-1i*k_x);
H(1, 4) = gamma_y+lambda_y*exp(1i*k_y);
H(2, 3) = (-gamma_y-lambda_y*exp(-1i*k_y))*(delta_x+1);
H(3, 1) = conj(H(1, 3));
H(4, 2) = conj(H(2, 4));
H(4, 1) = conj(H(1, 4));
H(3, 2) = conj(H(2, 3));
% 
BBH = Htrig(H);
BBH_sym = sym(BBH)
BBH = BBH.Subsall('sym');
BBH = BBH <'POSCAR_4';
BBH = BBH <'KPOINTS_4';
%%
delta_x = 0;
delta_y = 0;
Delta_x = 0;
Delta_y = 0;
delta = 0;
lambda = 1;
gamma = 0.5;
% bandplot
Fig = Figs(2,2);
BBH_n = BBH.Subsall();
bandplot(BBH_n.EIGENCAR_gen(),'title',"BBH-TETRA-Htrig",'KPOINTS','KPOINTS_4','ax',Fig.axes(1));
%%
[WccCAR,~,klist_l] = BBH_n.WannierCenter('script','nu_x(k_y)');
vasplib_plot.WccPlot(WccCAR,klist_l,'ax',Fig.axes(3),'title','\nu_x(k_y)','xlabel','k_y');
%% 
[nested_BFCAR,nested_klist_l] = BBH_n.nestedWilsonLoop('script','p_y(nu_x)','knum_evol',101,'knum_int',101,'V',diag([1,-1,-1,1]*1E-8));
 vasplib_plot.WccPlot(nested_BFCAR,nested_klist_l,'ax',Fig.axes(2),'xlabel','k_x');
title(Fig.axes(2),['$p_y^{\nu_x}= \frac{1}{N_{k_x}}\sum_{k_x} p_y^{\nu_x}(k_x)$ = ',num2str(mean(nested_BFCAR))],'Interpreter','latex');
%% 
[nested_BFCAR,nested_klist_l] = BBH_n.nestedWilsonLoop('script','p_x(nu_y)','knum_evol',101,'knum_int',101);
 vasplib_plot.WccPlot(nested_BFCAR,nested_klist_l,'ax',Fig.axes(4),'xlabel','k_y');
title(Fig.axes(4),['$p_x^{\nu_y}= \frac{1}{N_{k_y}}\sum_{k_y}p_x^{\nu_y}(k_x)$ = ',num2str(mean(nested_BFCAR))],'Interpreter','latex');
%% 
Fig.axes(1).FontSize = 16;
Fig.axes(2).FontSize = 16;
Fig.axes(3).FontSize = 16;
Fig.axes(4).FontSize = 16;