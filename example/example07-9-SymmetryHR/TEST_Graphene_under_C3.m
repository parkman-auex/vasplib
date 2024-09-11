clear
sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

Graphene = HR(2);
Graphene = Graphene<'POSCAR';
Graphene = Graphene<'KPOINTS';
search_range = [1 1 0];
maxR = 2.5;
Accuracy = 2;
Graphene = Graphene.nn(search_range, Accuracy ,maxR);
[Rnn,~,~,~] = Graphene.nn_information();
level_cut = 1;
Graphene = Graphene.H_TBSK_gen('level_cut',level_cut,'per_dir',[1 1 0]);
% list mode
Graphene = rewrite(Graphene);
Graphene.symvar_list;
% t = sym('t',[6,1],'real');
% Graphene.HcoeL = t;
list(Graphene);
Graphene.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
C3z = Oper.rotation(1/3,[0,0,1]);
C3z.U = eye(2);
%%
Graphene_C3 = Graphene.applyR(C3z.R);
list(Graphene_C3);
Graphene_C3.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
% [Graphene_test,VectorDistMat] = dualizeOper(Graphene,C3z);
%%
% Graphene_test = applyRU(Graphene,C3z );
% Graphene_test.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
C6z = Oper.rotation(1/6,[0,0,1]);
C6z.U = [0,1;1,0];
Graphene_C6 = Graphene.applyRU(C6z );
list(Graphene_C6);
Graphene_C6.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%