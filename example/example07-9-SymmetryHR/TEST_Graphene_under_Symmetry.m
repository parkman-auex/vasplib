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
Graphene = Graphene.nn_sk_smart(search_range, Accuracy ,maxR);
[Rnn,~,~,~] = Graphene.nn_information();
level_cut = 1;
Graphene = Graphene.H_TB_gen_SK(level_cut,[1 1 0]);
% list mode
Graphene = rewrite(Graphene);
Graphene.symvar_list;
t = sym('t',[6,1],'real');
Graphene.HcoeL = t+1i*t;
list(Graphene);
Graphene.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
Tr = Oper.time_reversal(3,eye(2));
%%
Graphene_Tr = Graphene.applyRU(Tr);
list(Graphene_Tr);
Graphene_Tr.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
syms E_s real;
Graphene_t = Graphene.set_hop_mat(sigma_0*(E_s+1i*E_s),[0,0,0],'sym');
list(Graphene_t);
Ch = Oper.chiral(3,sigma_z );
Graphene_Ch  = Graphene_t.applyRU(Ch);
list(Graphene_Ch);