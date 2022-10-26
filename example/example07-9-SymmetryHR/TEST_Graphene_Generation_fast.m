%% 
clear;
Graphene  = HR(2);
Graphene = Graphene < 'POSCAR' ;
%Graphene.Rm = double(Graphene.Rm);
%Graphene.orbL = double(Graphene.orbL);
Graphene= Graphene.nn([1,1,0],1e-4,1.15);
sym(Graphene.nn_store)
%%
Graphene2 = Graphene.init('level_cut',1,"onsite",1,'fast',true); % HR.init
%% better in live script
list(Graphene2);
% Graphene2.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
C3 = Oper.rotation(1/3,[0,0,1],0,eye(2));
Tr = Oper.time_reversal(3,diag([1 1]));
Chiral = Oper.chiral(3,diag([1,-1]));
Groups = generate_group([C3,Tr,Chiral]);
%%
Graphene_test = Graphene2.applyOper([Tr],'generator',true,'fast',true);
Graphene_test = Graphene_test.applyOper([C3],'generator',true,'fast',true);
Graphene_test = Graphene_test.applyOper([Chiral],'generator',true,'fast',true);

%%
Graphene_test = Graphene2.applyOper([C3,Tr,Chiral],'generator',true,'fast',true);
Graphene_test2 = Graphene_test.GenfromOrth();
list(Graphene_test2);
%Graphene_test2.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
Htrig_sym  = sym(Graphene_test2.HR2Htrig);
%%