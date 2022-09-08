clear
sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

Graphene = HR(2);
Graphene = Graphene<'POSCAR';
Graphene = Graphene<'KPOINTS';
%%
search_range = [1 1 0];
maxR = 2.5;
Accuracy = 2;
Graphene = Graphene.nn_sk_smart(search_range, Accuracy ,maxR);
[Rnn,~,~,~] = Graphene.nn_information();
level_cut = 1;
Graphene = Graphene.H_TB_gen_SK(level_cut,[1 1 0]);
%% list mode
Graphene_list = rewrite(Graphene);
Graphene_list.symvar_list
VppP_1 = 1;
list(Graphene_list)
Graphene_n = Graphene_list.Subsall();
EIGENCAR = Graphene_n.EIGENCAR_gen();
bandplot(EIGENCAR ,[-3,3]);