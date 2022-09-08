clear;
%% ThirdTerm 2D hex
%% 
% 
% useful tool

s_0   = pauli_matric(0)  ;  s_x = pauli_matric(1);  s_y =  pauli_matric(2) ;  s_z = pauli_matric(3);
sigma_0 = pauli_matric(0);sigma_x =  pauli_matric(1);sigma_y =  pauli_matric(2);sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0);  tau_x =  pauli_matric(1);  tau_y =  pauli_matric(2);  tau_z = pauli_matric(3);

syms B A0 A1 A2 real;
syms k_x k_y k_z real;
% Models
%
A       = A0+A2*(k_x^2+k_y^2)+A1*(k_z^2);
ThirdTerm = HK(4,3);
% 4: orbital number  dimension of mat
% 3: To k^3
%
ThirdTerm = ThirdTerm ...
    +Term(A*k_x ,sigma_y*tau_x )...
    +Term(-A*k_y ,sigma_x*tau_x )...    
    +Term(B*k_z*(k_x^2-k_y^2),sigma_x*tau_x )...
    +Term(-2*B*k_z*(k_x*k_y),sigma_y*tau_x )... 
;
ThirdTerm = ThirdTerm.Subsall('sym');
ThirdTerm = ThirdTerm <'POSCAR_Bi2Se3';
%Tr = Oper.time_reversal(3,double(-1i*gamma_matric(4,5)));% tau_z sigma_0
Mx = Oper.mirror([1,0,0],double( sigma_x*tau_0));%  
Tr = Oper.time_reversal(3,double(-1i*sigma_y*tau_0),nan);
% %My = Oper.mirror([0,1,0],double(1i*sigma_y*tau_z),nan);
I = Oper.inversion(3,diag([1,-1,1,-1]));
%C6z = Oper.rotation(1/6,[0,0,1],false,   double(expm(1i*sym(pi)*diag([1 3 -1 -3]/6)))); 

% C2z = Oper.rotation(1/2,[0,0,1],false,   expm(1i*(pi/2)*diag([1 3 -1 -3]))); 
C3z = Oper.rotation(1/3,[0,0,1],false,   double(expm(1i*sym(pi/3)*diag([-1 -1 1 1])))); 
% E = Oper.identity(3,4);
% %C4z = Oper.rotation(1/4,[0,0,1],false,sym(expm(1i*pi/4 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true)
% % C6z = Oper.rotation(1/6,[0,0,1],false,sym(expm(1i*pi/6 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true);
% % C3z = Oper.rotation(1/3,[0,0,1],false,sym(expm(1i*pi/3 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true);
% groups = [Tr,My,I,C6z];
%%
[ThirdTerm_Problem,ThirdTerm_trig]= ThirdTerm.kp2TB([0,0,0],[Mx,Tr,I,C3z],'level_cut',2,'Accuracy',1e-7,'mini',true);
ThirdTerm_TB = ThirdTerm_Problem.HR;
%ThirdTerm_trig.HcoeL = subs(ThirdTerm_trig.HcoeL,ThirdTerm_TB.symvar_list,[B/(-1.70279) sym(zeros(1,8))]);
%ThirdTerm_trig = ThirdTerm_trig.simplify();
%ThirdTerm_TB =subs(ThirdTerm_TB,ThirdTerm_TB.symvar_list,[B/(-1.70279) sym(zeros(1,8))]);
%ThirdTerm_TB = ThirdTerm_TB.simplify();
%%
ThirdTerm_TB.list()
vpa(ThirdTerm_TB.HcoeL)
vpa((ThirdTerm_trig.sym()))
ThirdTerm_trig_hk = ThirdTerm_trig.Htrig2HK('Order',2);
ThirdTerm.sym()
vpa(simplify(ThirdTerm_trig_hk.sym()))
