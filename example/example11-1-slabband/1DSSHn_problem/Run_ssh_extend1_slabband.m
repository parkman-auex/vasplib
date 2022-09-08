%% Run fot test cut_piece for ssh model
import vasplib_tool.*;
%% init
surfmode =1;
slabmode = 1;
repeatnum = 21 ;
norb_enforce = -1;
fermi_energy = 0;
glue_edges = 'false';
orbital_init = [0 0.5 0.5;0.4 0.5 0.5];
kpoint_single = [0 0 0];
%load system
if 1 ==1
    %!cp POSCAR_SSH POSCAR
    %!cp KPOINTS_1D KPOINTS

    POSCAR_read;

    Model_ssh_extend1;
    PARA_ssh_extend1;


    SSH_extend_n = Subsall(SSH_extend);
end
% bandplot
disp('bandplot');
tic;
    EIGENCAR = EIGENCAR_gen(SSH_extend_n,'t');
    bandplot(EIGENCAR,[-2,2]);
toc;

% surf_begin
H_test_for_surf = SSH_extend_n;
[temp_sparse,~,~] = hr2hr_sparse(H_test_for_surf);

% surf 
if surfmode ==1
    [temp_slab_sparse,fin_orb] = cut_piece(temp_sparse,repeatnum,1,glue_edges,orbital_init,'hrsn-sgl');
    if slabmode ==1
        [EIGENCAR_mono,WAVECAR_slab]  = EIGENCAR_gen(temp_slab_sparse,'ms',fermi_energy,norb_enforce,kpoint_single);
    end
end
% plot 

[fig,ax]=creat_figure();Y1 = EIGENCAR_mono;plot(ax,Y1,'-o');
orb_list = fin_orb;
[norb,~] = size(orb_list);
WaveFunc1 = WAVECAR_slab(:,norb/2);
WaveFunc2 = WAVECAR_slab(:,norb/2+1);
WaveFunc = [WaveFunc1 WaveFunc2 ];
 [~ , ~]=PARCHG_gen(orb_list,WaveFunc/8);

 % adjust parm
v_list = linspace(0.5,2,w_number);
EIGENCAR=[];
    for i=1:length(v_list)
        v = v_list(i);
        
        SSH_extend_n = Subsall(SSH_extend);
        hr_sparse = hr2hr_sparse(SSH_extend_n);
        [temp_slab_sparse,fin_orb] = cut_piece(hr_sparse,repeatnum,1,glue_edges,orbital_init,'hrsn-sgl');
        [EIGENCAR_mono,~]  = EIGENCAR_gen(temp_slab_sparse,'ms',fermi_energy,norb_enforce,kpoint_single);
        EIGENCAR=[EIGENCAR,EIGENCAR_mono];
    end
[Nbands,~]= size(EIGENCAR);    
[fig,ax]=creat_figure();
    for Ei=1:Nbands
        plot(ax,v_list,EIGENCAR(Ei,:),'LineWidth',1.0,'Color',[0 0 1],'DisplayName',num2str(Ei));
    end