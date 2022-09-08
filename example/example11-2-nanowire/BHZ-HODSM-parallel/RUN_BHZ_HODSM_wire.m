%% RUN BHZ_kzkplus model
% profile on -history;
clear;
% mode panel
bulkmode =2;
wiremode = 1;
nanodisk_mode =2;
%% load MODEL
!cp POSCAR_QWZ_2 POSCAR
MODEL_BHZ_test;
if 1==1 % 1 TB 2 kp
    K_xyz2K_123;
    Discretize_kp2TB;
    Split_Hamilton_kp2TB_smart;
end
%% para
PARA_QWZ;
B1 =   0;
% for wire band
norb_enforce = 1-1;
fermi = 0;

Nslab = [40,40,1];
np = 2; % mpi

WAN_NUM = 4;
orbital_init = ones(WAN_NUM,3)*0.5; % just the default set 
% you can use
% [EIGENCAR_wire,WEIGHTCAR_wire,orb_list]  =
% EIGENCAR_gen_wire(H_hr,fin_dir,np,fermi,norb_enforce)
%% subs
BHZ_HODSMn = simple_hr(Subsall(H_xyz,'struct'));
H_hr = hr2hr_sparse(BHZ_HODSMn);
%% plot
if bulkmode == 1
    EIGENCAR = EIGENCAR_gen(BHZ_HODSMn,'t');
    %bandplot(EIGENCAR,[-3,3],'BHZ_HODSM');
end

% [H_wire,Hnum_list_wire,vector_list_wire] = Hnanowire_gen(H_hr.HnumL,H_hr.vectorL,Nslab);
% [EIGENCAR_wire,~]  = EIGENCAR_gen(H_wire,'ts');
%  bandplot(EIGENCAR_wire,[-1,1],'Wire-BHZ-HODSM');
if wiremode == 1
    disp('direct wireplot');
    !cp KPOINTS_1D_Z_half KPOINTS
    tic
    %[EIGENCAR_wire,~]  = EIGENCAR_gen(temp_slab_sparse2,'ts');
    [EIGENCAR_wire,WEIGHTCAR_wire,orb_list]  = EIGENCAR_gen_wire(H_hr,Nslab,np,fermi,norb_enforce,orbital_init);
    toc
    
    tic
%     HSVCAR_surf = HSVCAR_gen(orb_list,'surf');
%     HSVCAR_hinge = HSVCAR_gen(orb_list,'hinge');
%     [COLORCAR_surf,WEIGHTCAR_surf] = COLORCAR_gen(WAVECAR_wire,HSVCAR_surf);
%     [COLORCAR_hinge,WEIGHTCAR_hinge] = COLORCAR_gen(WAVECAR_wire,HSVCAR_hinge);
    %bandplot(EIGENCAR_wire,[-1,1],'Wire-BHZ-HODSM');
    %WEIGHTCAR_wire_cell{1} = WEIGHTCAR_wire/100;
    EIGENCAR_wire = [fliplr(EIGENCAR_wire),EIGENCAR_wire];
    WEIGHTCAR_wire = [fliplr(WEIGHTCAR_wire),WEIGHTCAR_wire];
    !cp KPOINTS_1D_Z KPOINTS
    %pbandplot(WEIGHTCAR_wire,EIGENCAR_wire,[-1,1],'Wire-BHZ-HODSM-fat-hinge');
%     pbandplot(WEIGHTCAR_surf,EIGENCAR_wire,[-1,1],'Wire-BHZ-HODSM-fat-surf');
    mkdir('data_results');
    save('data_results/BHZ-HODSM-parallel','-v7.3');
    toc;
end
% profsave;
% profile off;