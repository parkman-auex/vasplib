clear;
%% Bi2Se3 2D hex
%% 
% 
% useful tool

s_0   = pauli_matric(0)  ;  s_x = pauli_matric(1);  s_y =  pauli_matric(2) ;  s_z = pauli_matric(3);
sigma_0 = pauli_matric(0);sigma_x =  pauli_matric(1);sigma_y =  pauli_matric(2);sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0);  tau_x =  pauli_matric(1);  tau_y =  pauli_matric(2);  tau_z = pauli_matric(3);

syms C0 C1 C2 real;
syms M0 M1 M2 real;
syms A0 A1 A2 B D real;
syms k_x k_y k_z real;
% Models

M       = M0-M2*(k_x^2+k_y^2)-M1*(k_z^2);
E0k     = C0+C2*(k_x^2+k_y^2)+C1*(k_z^2);
%A       = A0+A2*(k_x^2+k_y^2)+A1*(k_z^2);
A       = A0;
k_plus  = k_x + 1i* k_y;
k_minus = k_x - 1i* k_y;
%
Bi2Se3 = HK(4,2);
% 4: orbital number  dimension of mat
% 3: To k^3
%
Bi2Se3 = Bi2Se3 ...
    +Term(D*k_z,sigma_0*tau_y)...
    +Term(A*k_x ,sigma_y*tau_x )...
    +Term(-A*k_y ,sigma_x*tau_x )...
    +Term(E0k   ,sigma_0*tau_0 )...
    +Term(M     ,sigma_0*tau_z )...
;
Bi2Se3 = Bi2Se3.Subsall('sym');
Bi2Se3 = Bi2Se3 <'POSCAR_Bi2Se3';
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
[Bi2Se3_TB,Bi2Se3_trig]= Bi2Se3.kp2TB([0,0,0],[Mx,Tr,I,C3z],'level_cut',2,'Accuracy',1e-7);

%%
Bi2Se3_TB.list()
vpa(Bi2Se3_TB.HcoeL)
vpa((Bi2Se3_trig.sym()))
Bi2Se3_trig_hk = Bi2Se3_trig.Htrig2HK('Order',2);
Bi2Se3.sym()
vpa(simplify(Bi2Se3_trig_hk.sym()))
Bi2Se3_trig_hk_3 = Bi2Se3_trig.Htrig2HK('Order',3);
vpa(simplify(Bi2Se3_trig_hk_3.sym()))