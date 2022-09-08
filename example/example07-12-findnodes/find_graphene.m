%% initial
clear;
sigma_0 = [1 0;0 1]; sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0]; sigma_z = [1 0;0 -1];

Graphene = HR(2);
Graphene = Graphene<'POSCAR';
Graphene = Graphene<'KPOINTS';
search_range = [1 1 0];
maxR = 2.5;
Accuracy = 1e-2;
Graphene = Graphene.nn(search_range, Accuracy ,maxR);
Graphene = Graphene.H_TBSK_gen('level_cut',1,'per_dir',[1 1 0]);
%%
VppP_1 = 1;
Graphene_n = Graphene.Subsall();
% EIGENCAR = Graphene_n.EIGENCAR_gen();
% bandplot(EIGENCAR,[-3,3],'Color','blue');
%%
nk = [10,10,1];
[nodes_s, nodes_r] = findnodes(Graphene_n,'nk',nk);
%% 
Graphene_n.BZplot('mode','2D'); % 3D
nodes_r =  nodes_s*Graphene_n.Gk;
scatter3(nodes_r(:,1),nodes_r(:,2),nodes_r(:,3),100,'filled');