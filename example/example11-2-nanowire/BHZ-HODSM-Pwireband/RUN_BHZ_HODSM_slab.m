%% RUN BHZ_kzkplus model
clear;
% mode panel
bulkmode =2;
slabmode = 1;
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
glue_edges = 'false';
repeatnum = 10;
fin_direction = 2;
fin_direction2 =1;
norb_enforce = -1;
%% subs
BHZ_HODSM_2Dn = simple_hr(Subsall(H_xyz,'struct'));
Gen_hr(simple_hr(BHZ_HODSM_2Dn),'BHZ_HODSM_2Dn.dat');
%% plot
if bulkmode == 1
    EIGENCAR = EIGENCAR_gen(BHZ_HODSM_2Dn,'t');
    bandplot(EIGENCAR,[-3,3],'BHZ_HODSM-2D');
end
%% Slab
BHZ_HODSM_2Dn = simple_hr(BHZ_HODSM_2Dn);
[temp_sparse,~,~] = hr2hr_sparse(BHZ_HODSM_2Dn);
if slabmode == 1
 
    orbital_init =   [ 0.5 0.5 0.50000  ;...  % AgAs1 px up
                       0.5 0.5 0.50000  ;...  % AgAs2 px up
                       0.5 0.5 0.50000  ;...  % AgAs1 px dn]8907hjkjyui
                       0.5 0.5 0.50000  ];... % AgAs2 py dn   
    !cp POSCAR_QWZ_2 POSCAR
    tic
    disp('first direction');
    [temp_slab_sparse,fin_orb1] = cut_piece(temp_sparse,repeatnum,fin_direction,glue_edges,orbital_init,'hrsn-sgl');
    orb_list =  fin_orb1;
    toc
    !cp POSCAR_super_fin POSCAR
    tic
    !cp KPOINTS_1D_Z KPOINTS
%    [EIGENCAR_slab,WAVECAR_slab]  = EIGENCAR_gen(temp_slab_sparse,'ms');
%     HSVCAR_surf = HSVCAR_gen(orb_list,'surf');
%     HSVCAR_hinge = HSVCAR_gen(orb_list,'hinge');
%     [COLORCAR_surf,WEIGHTCAR_surf] = COLORCAR_gen(WAVECAR_slab,HSVCAR_surf);
%     [COLORCAR_hinge,WEIGHTCAR_hinge] = COLORCAR_gen(WAVECAR_slab,HSVCAR_hinge);
%     bandplot(EIGENCAR_slab,[-1,1],'Slab-BHZ-HODSM');
%     pbandplot(WEIGHTCAR_hinge,EIGENCAR_slab,[-1,1],'Slab-BHZ-HODSM-fat-hinge');
%     pbandplot(WEIGHTCAR_surf,EIGENCAR_slab,[-1,1],'Slab-BHZ-HODSM-fat-surf');
    % Gen_hr
    %BHZ_kzkplus_slab = hr_sparse2hr(temp_slab_sparse);
    %BHZ_kzkplus_slab = simple_hr(BHZ_kzkplus_slab);
    %[hrdat,splen] = Gen_hr(BHZ_kzkplus_slab,'wannier90_hr.dat.slab.sparse','sparse') ;
end

if wiremode == 1
    disp('second direction');
    [temp_slab_sparse2,fin_orb2] = cut_piece(temp_slab_sparse,repeatnum,fin_direction2,glue_edges,fin_orb1,'hrsn-sgl');
    !cp KPOINTS_1D_Z KPOINTS
    !cp POSCAR_super_fin POSCAR
    orb_list = fin_orb2;
    tic
    disp('EIGENCAR_caculate');
    HSVCAR = HSVCAR_gen(orb_list,'hinge');
    %
    save('HSVCAR.mat','HSVCAR');% mast
    %[EIGENCAR_wire,~]  = EIGENCAR_gen(temp_slab_sparse2,'ts');
    [EIGENCAR_wire,WEIGHTCAR_wire]  = EIGENCAR_gen(temp_slab_sparse2,'ts-wire');
    toc
    tic
%     HSVCAR_surf = HSVCAR_gen(orb_list,'surf');
%     HSVCAR_hinge = HSVCAR_gen(orb_list,'hinge');
%     [COLORCAR_surf,WEIGHTCAR_surf] = COLORCAR_gen(WAVECAR_wire,HSVCAR_surf);
%     [COLORCAR_hinge,WEIGHTCAR_hinge] = COLORCAR_gen(WAVECAR_wire,HSVCAR_hinge);
    bandplot(EIGENCAR_wire,[-1,1],'Wire-BHZ-HODSM');
    pbandplot(WEIGHTCAR_wire,EIGENCAR_wire,[-1,1],'Wire-BHZ-HODSM-fat-hinge');
%     pbandplot(WEIGHTCAR_surf,EIGENCAR_wire,[-1,1],'Wire-BHZ-HODSM-fat-surf');
    toc;
end





if nanodisk_mode ==1
    disp('second direction');

    [temp_slab_sparse2,fin_orb2] = cut_piece(temp_slab_sparse,repeatnum,fin_direction2,glue_edges,fin_orb1,'hrsn-sgl');
    toc
    !cp KPOINTS_1D KPOINTS
    !cp POSCAR_super_fin POSCAR
    orb_list = fin_orb2;
    
    % cut hexagonal disk
    % a pall nod
    tic
    disp('EIGENCAR_caculate');
    [EIGENCAR_mono,WAVECAR_disk]  = EIGENCAR_gen(temp_slab_sparse2,'ms',0,norb_enforce,[-0.00,0,0.0]);
    Y1 = EIGENCAR_mono;
    figure();
    % save('EuAgAs_kp1_soc1_EIGENCAR2_WAVECAR.mat','EIGENCAR2','WAVECAR_disk');
    plot(Y1,'-o');
    
    if norb_enforce <0
        [norb,~] = size(orb_list);
        xlim([norb/2-48 norb/2+48]);
    else
        norb = norb_enforce;
    end           
    WaveFunc1 = WAVECAR_disk(:,norb/2-2);
    WaveFunc2 = WAVECAR_disk(:,norb/2-1);
    WaveFunc3 = WAVECAR_disk(:,norb/2);
    WaveFunc4 = WAVECAR_disk(:,norb/2+1);
    WaveFunc5 = WAVECAR_disk(:,norb/2+2);
    WaveFunc6 = WAVECAR_disk(:,norb/2+3);
    WaveFunc = [WaveFunc3 WaveFunc4 WaveFunc1 WaveFunc2 WaveFunc5 WaveFunc6];
    % [~ , ~]=PARCHG_gen(orb_list,WaveFunc2);
    % [~ , ~]=PARCHG_gen(orb_list,WaveFunc3);
    % [~ , ~]=PARCHG_gen(orb_list,WaveFunc4);
    % [~ , ~]=PARCHG_gen(orb_list,WaveFunc5);
    [~ , ~]=PARCHG_gen(orb_list,WaveFunc);   
    toc
end

