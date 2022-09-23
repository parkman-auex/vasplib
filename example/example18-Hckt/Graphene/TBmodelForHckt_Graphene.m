%% Graphene
Type = 'list'; %listmat
syms v real;
Graphene = HR.from_POSCAR_SE('POSCAR','level_cut',1,'onsite',false,'per_dir',[1,1,0],'WAN_NUM',2,'r_max',3,'Type',Type);
Graphene = subs(Graphene,[v]);
Graphene = Graphene<'POSCAR';
%% prepare
Hcktpre = Graphene.HRforHckt();
% Add Unit
L_0 = 1e-6;
MU = 100;
try
    switch magnitude
        case 'p'
            Cfactor = MU*1E-12;% 100 pF
        case 'n'
            Cfactor = MU*1E-09;% 100 nF
        case 'u'
            Cfactor = MU*1E-06;% 100 uF
        case 'm'
            Cfactor = MU*1E-03;% 100 mF
    end
catch
    Cfactor = 100*1E-12;% 100 pF
    magnitude = 'p';
end
OmegaFactor = (Cfactor*L_0*100)^(-1/2);
%% Bulk:Armchair primitive
%
%HexTaTb_1 = HexTaTb_1.rewind();
% bulk
FigTest = Figs(1,2);
v= 1;
C_0 = 1;
%drawnow;
bandplot(Graphene.Subsall().EIGENCAR_gen('convention','I'),...
    [-3,3],'title',"t = "+num2str(v),'Color','b','ax',FigTest.axes(1));
%
Hcktpre_n = Hcktpre.Subsall();
EIGENCAR = Hcktpre_n.EIGENCAR_gen();
F0CAR = (EIGENCAR*L_0*Cfactor).^(-1/2)./(2*pi);
OMEGACUTpre = [0,2];
OMEGACUT = OMEGACUTpre*OmegaFactor;
bandplot(F0CAR,OMEGACUTpre*OmegaFactor,'ylabel','Frequency (Hz)','title',[...
    string(['C_v = ',num2str(v),magnitude,'F','; C_0 = ',num2str(C_0),magnitude,'F']);...
    ],...
    'ax',FigTest.axes(2));
%FigTest.Position = [];
return;
%HexTaTb_n.Gen_hr();
%!/Users/parkman/Documents/soft/wannier_tools-master/bin/wt.x

%% Slab
para_list = [1;2;3]
% kz = 0 armchair-2
for i = 1:size(para_list,1)
    v= para_list(i,1);
    HexTaTb_n = Graphene.Subsall();
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
    % plot
    bandplot(RCI_slab.EIGENCAR_gen(),[-3,3],klist_l,kpoints_l,kpoints_name,'title',"Edge:t_b = "+num2str(w)+", t_a = "+num2str(v),'Color','r');
end
