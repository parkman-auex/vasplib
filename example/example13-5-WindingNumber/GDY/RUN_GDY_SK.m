%Run
clear;
search_range = [1 1 1];
Accuracy = 3 ; %%

per_dir = [1 1 1] ;
% d_max = 6 ;% ai
% recommend_label_cut = sum(Rnn<d_max);

level_cut = 12;

% model
% graphene
if 1 == 1
[Rm,sites,Atom_name,Atom_num]=POSCAR_read('POSCAR','vasp');

[Atom_store,nn_store,Rnn]=nn_SK(Rm,sites,search_range, Accuracy);

%%
save('nn_store.mat');
else
    load nn_store.mat;
end


[H,H_Graphene,orb_list]=H_TB_gen_SK(level_cut,nn_store,Rm,per_dir );

% parm
E_onsite_1 = 0;
E_onsite_2 = 0;

% VppP_0  = 2.7;
% VppS_0  = -0.48;
%Rnn
% a0 = 1.23 ; delta0 = 0.45255; % VppP
% 
% d0 = 3.35; 
% %delta0=0;
% for i = 1:level_cut 
%  VppP_n = "VppP_"+string(i);
%   VppS_n = "VppS_"+string(i);
% tmp_string =  VppP_n+"  = VppP_0*exp(-(Rnn(i)-a0)/delta0) ";
% tmp_string2 =   VppS_n+" = VppS_0*exp(-(Rnn(i)-d0)/delta0)  ";
% eval(tmp_string);
% eval(tmp_string2);
% end
% give parm
VppP_1 = 3.34;
VppP_2 = 2.99;
VppP_3 = 2.82;
VppP_4 = 2.71*0.6;

VppP_5 = 0.0;
VppP_6 = 0;
VppP_7 = 0;
VppP_8 = 0;
VppP_9 = 0;
VppP_10 = 0;
VppP_11 = 0;
VppP_12 = 0;


VppS_1 = 0;
VppS_2 = 0;
VppS_3 = 0;
VppS_4 = 0;
VppS_5 = 0;
VppS_6 = 0;
VppS_7 = 0;
VppS_8 = 0;
VppS_9 = 0;
VppS_10 = 0;

VppS_11 = 0.463*1.4;
VppS_12 = 0.42*1.4;

Graphene_n = Subsall(H_Graphene);
% bulk band
EIGENCAR = EIGENCAR_gen(Graphene_n,'t');
% plot
bandplot(EIGENCAR,[-3,3],'Graphdiyne','b');

% Gen hr
Gen_hr(Graphene_n,'wannier90_hr.dat_06t1');  % generate hr.dat
