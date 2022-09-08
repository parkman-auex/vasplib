%% 
clear;
%;
BHZ_C6  = HR(4);
BHZ_C6  = 'POSCAR_BHZ_C6'>BHZ_C6 ;
Accuracy = 1e-4;
BHZ_C6= BHZ_C6.nn([1,1,1],Accuracy,2);
%%
BHZ_C6 = BHZ_C6.init('level_cut',2,"onsite",1); %最近邻
% BHZ_C6 = BHZ_C6.init('level_cut',2,"onsite",1); %次近邻
%% better in live script
% list(BHZ_C6);
% BHZ_C6.show('HOPPING','scale', 2.4560000896,'atomscale',0.1);
%%
C6 = Oper.rotation(1/6,[0,0,1],false,   expm(1i*(pi/6)*diag([1 3 -1 -3])));% 
Tr = Oper.time_reversal(3,double(-1i*gamma_matric(4,5)));% tau_z sigma_0
I = Oper.inversion(3,diag([1 -1 1 -1]));% -tau_x sigma_0
%Mx = Oper.mirror([1,0,0],double( 1i*gamma_matric(2,5)));% 
My = Oper.mirror([0,1,0],double( 1i*gamma_matric(1,2)));%    
Groups = generate_group([C6,Tr,I,My]);
%%
%BHZ_C6_test = BHZ_C6.applyOper(Groups)
BHZ_C6_test = BHZ_C6.applyOper([C6,Tr,I,My],'generator',true);
%BHZ_C6_test = BHZ_C6.applyOper(generate_group([C6,Tr]));
[BHZ_C6_test2,Sublist,Unique_term] = BHZ_C6_test.unique();
Varlist = BHZ_C6_test2.symvar_list;
%%
% syms t lambda_SO E_pz M_1 real;
% BHZ_C6_test = subs(BHZ_C6_test,Varlist,[t,E_pz,M_1,lambda_SO]);
% BHZ_C6_test = simplify(BHZ_C6_test);
% list(BHZ_C6_test);
BHZ_C6_test2.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
% %%
BHZ_C6_trig = BHZ_C6_test2.HR2Htrig();
BHZ_C6_trig = BHZ_C6_trig.simplify();
BHZ_C6_trig_sym = simplify(sym(BHZ_C6_trig));
%%
BHZ_C6_test_hk = BHZ_C6_trig.Htrig2HK('Order',2);
[BHZ_C6_test_hk,Sublist2,Unique_term2] =  unique(BHZ_C6_test_hk);
BHZ_C6_test_hk_sym = sym(BHZ_C6_test_hk);
%%
BHZ_C6_test_hk_3rd = BHZ_C6_trig.Htrig2HK('Order',3);
[BHZ_C6_test_hk_3rd,Sublist3,Unique_term3] =  unique(BHZ_C6_test_hk_3rd);
BHZ_C6_test_hk_sym_3rd = sym(BHZ_C6_test_hk_3rd);
%%
BHZ_C6_test_hk_1st = BHZ_C6_trig.Htrig2HK('Order',1);
[BHZ_C6_test_hk_1st,Sublist4,Unique_term4] =  unique(BHZ_C6_test_hk_1st);
BHZ_C6_test_hk_sym_1st = sym(BHZ_C6_test_hk_1st);