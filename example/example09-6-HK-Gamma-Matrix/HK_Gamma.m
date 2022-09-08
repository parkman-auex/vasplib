first_load = 1;
if first_load == 1
%% useful tool
s_0   = pauli_matric(0)  ;  s_x = pauli_matric(1);  s_y =  pauli_matric(2) ;  s_z = pauli_matric(3);
sigma_0 = pauli_matric(0);sigma_x =  pauli_matric(1);sigma_y =  pauli_matric(2);sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0);  tau_x =  pauli_matric(1);  tau_y =  pauli_matric(2);  tau_z = pauli_matric(3);
%%
Tr = Oper.time_reversal(3,double(-1i*gamma_matric(4,5)));% tau_z sigma_0
My = Oper.mirror([0,1,0],double( 1i*gamma_matric(1,2)));%  
%Tr = Oper.time_reversal(3,double(-1i*sigma_y*tau_0),nan);
%My = Oper.mirror([0,1,0],double(1i*sigma_y*tau_z),nan);
I = Oper.inversion(3,diag([1,-1,1,-1]));
C4z = Oper.rotation(1/4,[0,0,1],false,   expm(1i*(pi/4)*diag([1 3 -1 -3])),nan,'sym',true); 

%C4z = Oper.rotation(1/4,[0,0,1],false,sym(expm(1i*pi/4 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true)
C6z = Oper.rotation(1/6,[0,0,1],false,sym(expm(1i*pi/6 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true);
C3z = Oper.rotation(1/3,[0,0,1],false,sym(expm(1i*pi/3 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true);
groups = [Tr,My,I,C4z];
groups2 = [Tr,My,I,C6z];
groups3 = [Tr,My,I,C3z];


%%
BHZ_building = HK(4,1);
BHZ_building = BHZ_building <= 'POSCAR';
BHZ_building_C4_1 = BHZ_building.applyOper(groups,'generator',true);
[BHZ_building_C4_1st,Sublist1,Unique_term1] =  unique(BHZ_building_C4_1);
BHZ_building_C4_1st_sym = sym(BHZ_building_C4_1st);
%%
BHZ_building = HK(4,2);
BHZ_building = BHZ_building <= 'POSCAR';
BHZ_building_C4_2 = BHZ_building.applyOper(groups,'generator',true);
[BHZ_building_C4_2nd,Sublist2,Unique_term2] =  unique(BHZ_building_C4_2);
BHZ_building_C4_2nd_sym = sym(BHZ_building_C4_2nd);
%%
BHZ_building = HK(4,3);
BHZ_building = BHZ_building <= 'POSCAR';
BHZ_building_C4_3 = BHZ_building.applyOper(groups,'generator',true);
[BHZ_building_C4_3rd,Sublist3,Unique_term3] =  unique(BHZ_building_C4_3);
BHZ_building_C4_3rd_sym = sym(BHZ_building_C4_3rd);
% save
end
%%
[BHZ_building_C4_1st_Gamma,BHZ_building_C4_1st_latex_Gamma] = BHZ_building_C4_1st.GammaDecomposition();
[BHZ_building_C4_1st_pauli,BHZ_building_C4_1st_latex_pauli] = BHZ_building_C4_1st.pauliDecomposition();
[BHZ_building_C4_2nd_Gamma,BHZ_building_C4_2nd_latex_Gamma] = BHZ_building_C4_2nd.GammaDecomposition();
[BHZ_building_C4_2nd_pauli,BHZ_building_C4_2nd_latex_pauli] = BHZ_building_C4_2nd.pauliDecomposition();
[BHZ_building_C4_3rd_Gamma,BHZ_building_C4_3rd_latex_Gamma] = BHZ_building_C4_3rd.GammaDecomposition();
[BHZ_building_C4_3rd_pauli,BHZ_building_C4_3rd_latex_pauli] = BHZ_building_C4_3rd.pauliDecomposition();