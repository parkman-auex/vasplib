%% RUN_HR_TBSK to generate 2D numeric wannier90_hr.dat  and maybe band
clear;
%% set
%Subs_mode  = 1;% numeric TB
Type = 'sparse';
% first_load = 1;
bulk_band  = 1;
slab_band  = 1;
conner_state  = 0;
% para
search_range = [1 1 0];
Accuracy = 3;
d_max_search = 4.4729;%4.4728;%2.410;% 1.4255;
d_max = d_max_search;% 1.4255;
Num = 1;  % max_p
N_rep = 2;
p =1;
%% para pre
d0 = 1.4249;
a0 = 1.424;
VppS_0  = -0.48;
VppP_0  = 2.727;%*10;%2.7;
delta0 = 0.45255;%0.39337;%0.45255;
para_struct.strvar = ["VppS","VppP"];
para_struct.numvar = [VppS_0,VppP_0];
para_struct.delta  = delta0 ;
%% H_TB gen for all Num
% [~,sites,~,~,~]=vasplib.POSCAR_readin('POSCAR_tb_1','tbsk');
% H_hr = HR(length(sites),'Type','sparse');
% H_hr = H_hr < 'POSCAR_tb_1';
% H_hr = H_hr.nn(search_range,Accuracy,d_max);
% H_hr = H_hr.H_TBSK_gen_sparse(...
%     'chiral',false,...
%     'onsite',false,...
%     'level_cut',3,...
%     'Rd',[d0,a0],...
%     'para',para_struct,...
%     'deltarule',1 ...
%     );
%%
%H_TBSK = HR.from_POSCAR_SE('POSCAR_1.08_SK',...
H_TBSK = HR.from_POSCAR_SE('POSCAR_tb_1',...
    'Type',Type,...
    'r_max',d_max,...
    'level_cut',3,...
    'search_range',search_range,...
    'Accuracy',1e-3,...
    'chiral',false,...
    'onsite',false,...
    'deltarule',1,...
    'Rd',[d0,a0],...
    'para',para_struct...
    );
%% band
SYSTEM = "band";
if bulk_band == 1
    % Efermi = 1.2;
    H_TBSK_n = H_TBSK < 'KPOINTS';
    EIGENCAR = H_TBSK_n.EIGENCAR_gen();
    % plot
    bandplot(EIGENCAR,[-1,1],'title',SYSTEM,'Color','b','POSCAR','POSCAR_tb_1');
end
%     %% Gen hr
%     % H_TBSK_n.Gen_hr('wannier90_hr.dat');  % generate hr.dat used for ir2tb
%% slab_band
tic;
if slab_band == 1 
    repeatnum   = 5; % too small?
    fin_dir     = 2; 
    glue_edges  = false;
    vacuum_mode = false;
    % H_TBSK_n = full(H_TBSK_n );
    % Gen Slab
    H_TBSK_n_slab = H_TBSK_n.cut_piece(repeatnum,fin_dir,glue_edges,vacuum_mode);
    % load KPOINTS
    H_TBSK_n_slab = H_TBSK_n_slab < 'KPOINTS_slab_X';
    % slab band
    Efermi = 0.0;
    EIGENCAR_slab = H_TBSK_n_slab.EIGENCAR_gen('LWAVE',false,'norb',-1);
    % plot
    bandplot(EIGENCAR_slab,[-1,1],'POSCAR','POSCAR_tb_1','KPOINTS','KPOINTS_slab_X');
end
toc;
%% Green func
fin_dir = 2;
KPOINTS_surf = 'KPOINTS_slab_X';
principle_layer = 1;
w_range = [-0.5,0.5,50];
eta =0.01;
[DOSCAR_l,~,~,w_list,klist_l,kpoints_l,kpoints_name] = ...
  H_TBSK_n.surf(w_range,fin_dir,KPOINTS_surf,principle_layer,eta);
heatplot(DOSCAR_l,w_list,klist_l,kpoints_l,kpoints_name);
%% corner state
if conner_state  == 1
    Nslab = [N_rep N_rep 1];
    np = 1;
    vacuum_mode = 1;
    H_TBSK_n_disk = H_TBSK_n.Hnanowire_gen(Nslab,np,vacuum_mode);
    % disk eigen func
    Efermi = 0;
    norb_enforce = -1;
    kpoint_s = [0 0 0];
    [EIGENCAR_disk,WAVECAR_disk,WEIGHTCAR] = H_TBSK_n_disk.EIGENCAR_gen('fermi',Efermi,'norb',norb_enforce,'klist',kpoint_s,'WEIGHTCAR',true);
    figure();
    norb = length(EIGENCAR_disk);
    fig_energy=scatter(1:norb,EIGENCAR_disk,1,...
        (WEIGHTCAR),'filled');
    %saveas(fig_energy,"energy_level_"+p);
    WaveFunc = WAVECAR_disk(:,4*(3*p^2+3*p+1)*N_rep*N_rep/2);
    orb_list = H_TBSK_n_disk.orbL;
    PARCHG_gen(orb_list,WaveFunc);
    view(0,90);
end






