%% 
clear;
%;
KaneMele  = HR(4);
KaneMele  = KaneMele<'POSCAR_KM' ;
KaneMele= KaneMele.nn([1,1,0],1e-2,1.15);
%%
KaneMele = KaneMele.init('level_cut',2,"onsite",1,'fast',true);
%% better in live script
% list(KaneMele);
% KaneMele.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
C3 = Oper.rotation(1/3,[0,0,1],false,double(expm(-1i*(pi/3)*gamma_matric(2,4))));% tau_0 e^(ipi/3*sigma_z)
Tr = Oper.time_reversal(3,double(-1i*gamma_matric(4,5)));% tau_0 i*sigma_y
I = Oper.inversion(3,double(-gamma_matric(1)));% -tau_x sigma_0
Mx = Oper.mirror([1,0,0],double( 1i*gamma_matric(2,5)));% not change sub lattice; tau_0 i*sigma_x
My = Oper.mirror([0,1,0],double( 1i*gamma_matric(2,3)));%     change sub lattice; tau_x i*sigma_y
Groups = generate_group([C3,Tr,I,Mx,My]);
%%
KaneMele_test = KaneMele.applyOper([C3,I,Mx,My,Tr],'generator',true,'fast',true);
KaneMele_test2 = KaneMele_test.GenfromOrth();
list(KaneMele_test2);
%Graphene_test2.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
Varlist = KaneMele_test2.symvar_list;
%%
syms t lambda_SO E_pz M_1 real;
KaneMele_test2 = subs(KaneMele_test2,Varlist,[lambda_SO,E_pz,t]);
KaneMele_test2 = simplify(KaneMele_test2);
KaneMele_test2.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
KaneMele_trig = KaneMele_test2.HR2Htrig();
KaneMele_trig_sym = simplify(sym(KaneMele_trig));