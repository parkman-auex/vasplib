%% 
clear;
Graphene  = HR(2);
Graphene = 'POSCAR'>Graphene;
Graphene.Rm
Graphene.orbL
Graphene= Graphene.nn([1,1,0],4,1.15);
sym(Graphene.nn_store)
%%
Graphene = Graphene.init('level_cut',1,"onsite",1);

%% better in live script
list(Graphene);
Graphene.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
C3 = Oper.rotation(1/3,[0,0,1],false,eye(2));
Tr = Oper.time_reversal(3,diag([1 1]));
Chiral = Oper.chiral(3,diag([1,-1]));
Groups = generate_group([C3,Tr,Chiral]);
%%
Graphene_test = Graphene.applyOper([C3,Tr,Chiral],'generator',true);
list(Graphene_test);
Graphene_test.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
[Graphene_test2,Sublist,Unique_term] = Graphene_test.unique();
Htrig_sym  = sym(Graphene_test2.HR2Htrig);
%%

Graphene_test3 = Graphene.applyOper(C3);
Graphene_test3 = Graphene_test3.applyOper(C3);
Graphene = Graphene_test3.applyOper(Tr);
% Graphene_test.Htrig_sym;