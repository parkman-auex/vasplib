%% HRclass Tutorial
% 基于HRclass实现石墨烯Kane-Mele model并 supercell (换边)
% 
% 2021.05
%% 
% * Author: 杨帆曾旭涛
% * Email: parkman@buaa.edu.cn
%% 

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
Accuracy = 1e-6;
Graphene = Graphene.nn(search_range, Accuracy ,maxR);
[Rnn,~,~,~] = Graphene.nn_information();
%% 
% 

Graphene = Graphene.H_TBSK_gen('level_cut',1,'per_dir',[1 1 0]);
Graphene.symvar_list
VppP_1 = 1;
list(Graphene);
%% Supercell_hr
% 通过一个变换矩阵可以对 一个HR对象进行 supercell（旋转+扩胞）

Ns = [1,0,0;1 2 0;0 0 1];
Graphene_semiarmchair_n = Graphene.Subsall.rewrite.supercell_hr(Ns);
Graphene_semiarmchair_n = Graphene_semiarmchair_n < 'KPOINTS_rectangle';
EIGENCAR = Graphene_semiarmchair_n.EIGENCAR_gen();
bandplot(EIGENCAR ,[-3,3],'title','Graphene-semiArmchair','Color','b','KPOINTS','KPOINTS_rectangle');
Graphene_semiarmchair_n.show('TwoD',true);axis equal
% Ns = [1,0,0;1 2 0;0 0 1];
Graphene_semiarmchair_n = [Graphene_semiarmchair_n;Graphene_semiarmchair_n];
Graphene_semiarmchair_n.quantumL(:,end) = [0.5;0.5;.5;0.5;-0.5;-0.5;-0.5;-0.5];
RANDPERM_GO = randperm(Graphene_semiarmchair_n.WAN_NUM);
%RANDPERM_GO = 1:8
Graphene_semiarmchair_n = Graphene_semiarmchair_n.reseq(RANDPERM_GO);
%% fold hr
%Graphene_fold_n = Graphene_semiarmchair_n.unfold_hr(Ns,"Accuracy",1e-4);
ORB_IDL = [1,2,1,2,3,4,3,4];
Graphene_unfold_n = Graphene_semiarmchair_n.unfold_hr(Ns,"Accuracy",1e-4,'orb_idL',ORB_IDL(RANDPERM_GO));
Graphene_unfold_n = Graphene_unfold_n.input_Rm('POSCAR');
Graphene_unfold_n = Graphene_unfold_n < 'KPOINTS';
EIGENCAR = Graphene_unfold_n.EIGENCAR_gen();
bandplot(EIGENCAR ,[-3,3],'title','Graphene-unfold','Color','b');
Graphene_unfold_n.show('TwoD',true);axis equal
%% Kane-Model model

V_NNN = [1 0 0; 0 1 0; -1 -1 0];
VppP_2 = 0.1 * VppP_1;

Kane_Mele = [Graphene,Graphene];
Kane_Mele = Kane_Mele.set_hop(...
     kron(sigma_z, -1i*VppP_2),1,1,V_NNN,'sym');
Kane_Mele = Kane_Mele.set_hop(...
     kron(sigma_z,  1i*VppP_2),2,2,V_NNN,'sym');
 
Kane_Mele = Kane_Mele.autohermi()
Kane_Mele_n = Kane_Mele.Subsall();
EIGENCAR = Kane_Mele_n.EIGENCAR_gen();
bandplot(EIGENCAR ,[-3,3]);
%% Supercell_hr
% 通过一个变换矩阵可以对 一个HR对象进行 supercell（旋转+扩胞）

Ns = [1,0,0;1 2 0;0 0 1];
Kane_Mele_semiarmchair_n = Kane_Mele_n.rewrite.supercell_hr(Ns);
Kane_Mele_semiarmchair_n = Kane_Mele_semiarmchair_n < 'KPOINTS_rectangle';
EIGENCAR = Kane_Mele_semiarmchair_n.EIGENCAR_gen();
bandplot(EIGENCAR ,[-3,3],'title','KaneMele-semiArmchair','Color','b','KPOINTS','KPOINTS_rectangle');
Kane_Mele_semiarmchair_n.show();
%%
Ns = [1,-1,0;1 2 0;0 0 1];
Kane_Mele_armchair_n = Kane_Mele_n.rewrite.supercell_hr(Ns);
Kane_Mele_armchair_n = Kane_Mele_armchair_n < 'KPOINTS';
EIGENCAR = Kane_Mele_armchair_n.EIGENCAR_gen();
% [klist_l,kpoints_l,kpoints_name] = Kane_Mele_armchair_n.kpath_information();
bandplot(EIGENCAR ,[-3,3],'title','KaneMele-Armchair','Color','r');
Kane_Mele_armchair_n.show();
%% Slab边缘态（surface）

repeatnum   = 20;
fin_dir     =  2;
glue_edges  = false;
vacuum_mode = 1;
% Gen Slab
Kane_Mele_n_slab = Kane_Mele_armchair_n.cut_piece(repeatnum,fin_dir,glue_edges,vacuum_mode);
% load KPOINTS
Kane_Mele_n_slab = Kane_Mele_n_slab < 'KPOINTS_slab';
% slab band
Efermi = 0.0;
% compare speed
tic;
EIGENCAR_slab = Kane_Mele_n_slab.EIGENCAR_gen();
toc;
Kane_Mele_armchair_n_slab_list = Kane_Mele_n_slab.rewrite();
tic
EIGENCAR_slab = Kane_Mele_armchair_n_slab_list.EIGENCAR_gen();
toc;
% plot
bandplot(EIGENCAR_slab,[-3,3],...
'KPOINTS','KPOINTS_slab',...
'title','KaneMele-Armchair-slab', ...
'Color','g');
%% Slab边缘态（using supercell_hr with OBC）
% Gen Slab
Kane_Mele_n_slab = Kane_Mele_n.supercell_hr(diag([1 20 1]),'Accuracy',1e-5,'OBC',[0 1 0]);
% load KPOINTS
Kane_Mele_n_slab = Kane_Mele_n_slab < 'KPOINTS_slab';
EIGENCAR_slab = Kane_Mele_n_slab.EIGENCAR_gen();
toc;
% plot
bandplot(EIGENCAR_slab,[-3,3],...
'KPOINTS','KPOINTS_slab',...
'title','KaneMele-slab(By supercell\_hr)', ...
'Color','r');

%% 能态与角落态
% 最后，我们将展现Kane-Mele的Armchair 边缘态；并使用 norb_enforce使用eigs加速

Nslab = [10 10 1];
np = 1;
vacuum_mode = 1;
% Gen disk-
n_disk = Kane_Mele_armchair_n.Hnanowire_gen(Nslab,np,vacuum_mode);
norb_enforce = 10;
fermi = 0;
[EIGENCAR_disk,WAVECAR_disk] = n_disk.EIGENCAR_gen('klist',[0 0 0],'LWAVE',true);
figure();
plot(EIGENCAR_disk,'-o');
ylim([-1 1]);
%% 
% armchair 边的边缘态有一点点问题?
WaveFunc = WAVECAR_disk(:,length(EIGENCAR_disk)/2:length(EIGENCAR_disk)/2+1);
orb_list = n_disk.orbL;
vasplib.waveplot(orb_list,WaveFunc);
axis equal
view(0,90)
%% 参考文献
% DOI: 10.1103/PhysRevLett.124.166804 
% 
%