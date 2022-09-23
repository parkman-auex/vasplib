%% SOTI naïve model
%% 
% 
% 
% 2D HexTaTb_Hcktpre model ta tb effect parkman 2022-01-12

% ------------------------    POSCAR    ------------------------
% TaTbHex
% 1.0
% 4.0000000000         0.0000000000         0.0000000000
% -2.0000000000         3.4641016151         0.0000000000
% 0.0000000000         0.0000000000         3.0000000000
% C Si
% 3  3
% Direct
% 0.300000012         0.000000000         0.000000000 C s 
% 0.000000000         0.300000012         0.000000000 C s 
% 0.699999988         0.699999988         0.000000000 C s 
% 0.699999988         0.000000000         0.000000000 Si s 
% 0.000000000         0.699999988         0.000000000 Si s 
% 0.300000012         0.300000012         0.000000000 Si s 
% ------------------------    ______    ------------------------
% %% Constuct TB model
% Type = 'list';%listmat
% HexTaTb = HR.from_POSCAR_SE('POSCAR','level_cut',2,'onsite',false,'per_dir',[1,1,0],'WAN_NUM',6,'r_max',3,'Type',Type);
% % 
% syms t_a t_b real;
% syms E1 E2 v w real;
% HexTaTb = HexTaTb.subs([E1 E2 t_a t_b]);
% HexTaTb_1 = subs(HexTaTb,[E1 E2 t_a t_b],[E1 -E1 v w]);
% HexTaTb_1 = HexTaTb_1<'POSCAR_sym';
%% armchair 2
Type = 'list';%listmat
syms v w real;
HexTaTb_2 = HR.from_POSCAR_SE('POSCAR_2','level_cut',2,'onsite',false,'per_dir',[1,1,0],'WAN_NUM',6,'r_max',3,'Type',Type);
HexTaTb_2 = subs(HexTaTb_2,[v w]);
HexTaTb_2 = HexTaTb_2<'POSCAR_sym_2';
%% prepare
HexTaTb_Hcktpre = HexTaTb_2;
HexTaTb_Hcktpre.HcoeL = -HexTaTb_Hcktpre.HcoeL;
% make symbolic term > 0
syms v w C_0 real;
assume(v>0);
assume(w>0);
assume(C_0>0);
BaseOnsiteL = repmat(C_0,[HexTaTb_Hcktpre.WAN_NUM,1]);
maxOnsite = C_0;
% set onsite
for i = 1:HexTaTb_Hcktpre.WAN_NUM
    maxOnsite = max(BaseOnsiteL(i)-sum(HexTaTb_Hcktpre.HcoeL(HexTaTb_Hcktpre.vectorL(:,4)==i)),maxOnsite);
end
for i = 1:HexTaTb_Hcktpre.WAN_NUM
    BaseOnsiteL(i) = maxOnsite + sum(HexTaTb_Hcktpre.HcoeL(HexTaTb_Hcktpre.vectorL(:,4)==i));
    HexTaTb_Hcktpre = HexTaTb_Hcktpre.set_hop(maxOnsite,...
            i,i,[0,0,0],'sym');
end
% Add Unit
L_0 = 1e-6;
MU = 100;
Cfactor = 100*1E-12;% 100 pF
OmegaFactor = (Cfactor*L_0*100)^(-1/2);
%% Bulk:Armchair primitive
%
%HexTaTb_1 = HexTaTb_1.rewind();
% bulk
para_list = [1,2];
for i = 1:size(para_list,1)
    FigTest = Figs(1,2);
    v= para_list(i,1);
    w= para_list(i,2);
    C_0 = 1;
    %drawnow;
    bandplot(HexTaTb_2.Subsall().EIGENCAR_gen('convention','I'),...
        [-3,3],'title',"t_b = "+num2str(w)+", t_a = "+num2str(v),'Color','b','ax',FigTest.axes(1));
    %
    HexTaTb_Hcktpre_n = HexTaTb_Hcktpre.Subsall();
    EIGENCAR = HexTaTb_Hcktpre_n.EIGENCAR_gen();
    F0CAR = (EIGENCAR*L_0*Cfactor).^(-1/2)./(2*pi);
    OMEGACUTpre = [0,2];
    OMEGACUT = OMEGACUTpre*OmegaFactor;
    bandplot(F0CAR,OMEGACUTpre*OmegaFactor,'ylabel','Frequency (Hz)','title',[...
        string(['C_A = ',num2str(MU*double(subs(BaseOnsiteL(1)))),'pF','; C_B = ',num2str(MU*double(subs(BaseOnsiteL(2)))),'pF']);...
        string(['C_v = ',num2str(MU*v),'pF','; C_w = ',num2str(MU*w),'pF']);...
        ],...
        'ax',FigTest.axes(2));
end

%HexTaTb_n.Gen_hr();
%!/Users/parkman/Documents/soft/wannier_tools-master/bin/wt.x

%% Slab

% kz = 0 armchair-2
for i = 1:size(para_list,1)
    v= para_list(i,1);
    w= para_list(i,2);
    HexTaTb_n = HexTaTb_Hcktpre_n.Subsall();
    repeatnum   = 30;
    fin_dir     =  2;
    glue_edges  = false;
    vacuum_mode = 0;
    % Gen Slab
    RCI_slab = HexTaTb_n.cut_piece(repeatnum,fin_dir,glue_edges,vacuum_mode);
    % load KPOINTS
    RCI_slab = RCI_slab < 'KPOINTS_slab';
    % slab band
    [klist_l,kpoints_l,kpoints_name] = RCI_slab.kpath_information();
    EIGENCAR = RCI_slab.EIGENCAR_gen();
    % plot
    F0CAR = (EIGENCAR*L_0*Cfactor).^(-1/2)./(2*pi);
    OMEGACUTpre = [0,2];
    OMEGACUT = OMEGACUTpre*OmegaFactor;
    bandplot(F0CAR,OMEGACUTpre*OmegaFactor,klist_l,kpoints_l,kpoints_name,'ylabel','Frequency (Hz)','title',"Edge:C_v = "+num2str(MU*v)+"pF, C_w = "+num2str(MU*w)+"pF",'Color','r');
end
%% Corner
Noccu = 3;
Nsuper = 30;
Nslab = [Nsuper Nsuper 1];
Nocc = Nslab(1)* Nslab(1) *Noccu*3/4;
np = 1;
vacuum_mode = 0;
% kz = 0 armchair-1
for i = 1:size(para_list,1)
    v = para_list(i,1);
    w = para_list(i,2);
    n_disk = HexTaTb_n.rewind.Hnanowire_gen(Nslab,np,vacuum_mode);
    n_disk_cut = n_disk;
    %         efsmall = 1e-6;
    %         rmfunc = @(i1,i2,i3) ...
    %             i1 > 1-1/Nsuper-2*TRACK+efsmall |...
    %             i2 > 1-1/Nsuper-2*TRACK+efsmall |...
    %             i2-i1 > 0.5- 0.5/Nsuper-TRACK+efsmall | ...
    %             i2-i1 < -0.5+0.5/Nsuper+TRACK-efsmall;
    %         n_disk_cut = n_disk.cut_orb('rmfunc',rmfunc);
    %n_disk_cut.show('NRPTS','vectorList',[0,0,0],'scale',5);
    %view(0,90)
    orb_list = n_disk_cut.orbL;
    norb_enforce = -1;
    fermi = 0;
    PM = 'select-points';
    minx = 0;miny = 0;
    maxx = 1;maxy = 1;
    %Center_plus = [minx miny 0.5;maxx maxy 0.5;];
    %Center_minus = [minx maxy 0.5;maxx miny 0.5];
    Center_plus = [minx miny 0.;maxx maxy/2 0.;maxx/2 maxy 0.;];
    Center_minus = [maxx maxy 0.;minx maxy/2 0.;maxx/2 miny 0.;];
    PS = struct('discrimination',0.1,'center',[Center_plus;Center_minus],'orientation',0,'sign',false);
    %error('as');
    [EIGENCAR_disk,WAVECAR_disk,WEIGHTCAR] = n_disk_cut.EIGENCAR_gen('WEIGHTCAR',true,'klist',[0 0 0]);
    %%
    figure = vasplib_tool.creat_figure();
    Nocc = length(EIGENCAR_disk)/2;
    %plot(EIGENCAR_disk,'-o');
    F0CAR_disk = (EIGENCAR_disk*L_0*Cfactor).^(-1/2)./(2*pi);
    scatter(1:length(EIGENCAR_disk),F0CAR_disk,100*ones(length(EIGENCAR_disk),1),(WEIGHTCAR).^3,'filled');
    colormap(ColorMap.Matplotlib('ocean'))
    xlim([Nocc-Nsuper Nocc+Nsuper]);
    title("Disk:C_v = "+num2str(MU*v)+"pF, C_w = "+num2str(MU*w)+"pF");
    %%
    WaveFunc = WAVECAR_disk(:,[Nocc:Nocc+1]);
    orb_list = n_disk_cut.orbL;
    vasplib.waveplot(orb_list,WaveFunc);
    axis equal
    view(0,90)

end

