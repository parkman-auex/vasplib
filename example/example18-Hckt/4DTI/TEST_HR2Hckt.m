%% Declare the primitive of 4D_TI
TITLE = '4D_TI';
Dim = 4;
%meshs = [3,3,3,3];
meshs = [60,2,3,3];
magnitude = 'p';
workdir = pwd;
Analysis = 'ac';
switch magnitude
    case 'p'
        OmegaCut = [0.4e7 0.8e7];
    case 'n'
        OmegaCut = [0 6e5];
    case 'u'
        OmegaCut = [0 2e4];
    otherwise
        OmegaCut = [0 2e7];
end
%%  run TB part
TBmodelForHckt_4D;
%% TRANSLATE
HCKT_test = Hcktpre_n.HR2Hckt('title',TITLE,'dim',Dim,'autominus',true,'magnitude',magnitude);
%%
HCKT_test.AC = '.AC LIN 400 0.4e7 0.8e7';
if 1 == 1
    %% Gen_sp
    Hckt.Genlib('magnitude',magnitude);
    [~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'analysis',Analysis,'ExciteNodes',[1;]);
    if strcmp(Analysis,'ac')
        DATAname = [char(Basename),'_AC','.ac0'];
        SPfile = strcat(Basename,'_AC.sp');
        EvalStringDATA = (['!hspiceAC ',DATAname]);
        SPLIS = strcat(Basename,'_AC.lis');
    else
        DATAname = [char(Basename),'.tr0'];
        SPfile = strcat(Basename,'.sp');
        EvalStringDATA = (['!hspiceTR ',DATAname]);
        SPLIS = strcat(Basename,'.lis');
    end
end
if isunix && ~ismac()
%% run; run it on linux system
EvalStringRun = (['!hspice -mt 20 -i ',SPfile,' |tee log']);
eval(EvalStringRun);
eval(EvalStringDATA);
elseif ispc()
    EvalStringRun = (['hspice -mt 4 -i ',SPfile,' -o ', SPLIS]);
    system(EvalStringRun);
    win_matlab.hspiceACTR(DATAname);
end
%% end
%% 
if (isunix && ~ismac) ||(ispc)
%% data clean
simulation_result = Hckt.read_hspice_ac(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result,'analysis',Analysis);
%% grid
if 1 == 1
SelectL = VectorList(:,end)==2;
ObservationsMat2 = ObservationsMat(SelectL,:);
VectorList2 = VectorList(SelectL,:);
else
    VectorList2 = VectorList;
    VectorList2(:,end) = 2;
    %ObservationsMat2 = ObservationsMat(SelectL,:);
    [VectorList2,ObservationsMat2] = HollowKnight.generalcontractrow2(VectorList2,ObservationsMat);
end
%% 
DOSCAR = Hckt.CollectVstruct(VectorList2,ObservationsMat2,OmegaCut,SpectrumX);
% [klist_cart,klist_frac,klist_l,kpoints_l,~] = ...
%     vasplib.kpathgen(...
%     [ ...
%      0 0 0 0; ...
%      1 0 0 0; ... 
%      0 0 1/3 -1/3; ...
%      0 1 1/3 -1/3; ...
%      0 5/12 0 -1/3; ...
%      0 5/12 1 -1/3; ...
%      0 5/12 1/3 0; ...
%      0 5/12 1/3 1; ...
%      ],...
%     60,(2*pi*eye(4)/eye(4)).','Dim',4);
% kpoints_name = ["\Gamma","\Gamma_x|\Gamma","\Gamma_y|\Gamma","\Gamma_z|\Gamma","\Gamma_w"];
[klist_cart,klist_frac,klist_l,kpoints_l,~] = ...
    vasplib.kpathgen(...
    [ ...
     0 0 0 0; ...
     1 0 0 0; ... 
     ],...
    60,(2*pi*eye(4)/eye(4)).','Dim',4);
kpoints_name = ["0","2\pi"];
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR,"klist_frac",klist_frac,"method","nearest");
Hckt.bandplot((EIGENCAR),OmegaCut,SpectrumX,...
    'title',['4D TI',' n = ',num2str(meshs)],...
    'klist_l',klist_l,'kpoints_l',kpoints_l,'kpoints_name',kpoints_name);
% savefig(gcf,'Fig3.fig');
% delete(gcf);

if 1 == 2
    %%
    SelectL = VectorList(:,end)==2;
    ObservationsMat2 = ObservationsMat(SelectL,:);
    VectorList2 = VectorList(SelectL,:);
    DOSCAR = Hckt.CollectVstruct(VectorList2,ObservationsMat2,OmegaCut,SpectrumX);
    EIGENCAR = zeros(400,60);
    for i = 1:numel(DOSCAR)
        EIGENCAR(i,:) = DOSCAR{i}(1,:,2,3).';
    end
    [klist_cart,klist_frac,klist_l,kpoints_l,~] = ...
        vasplib.kpathgen(...
        [ ...
        0 0 0 0; ...
        0 1 0 0; ...
        ],...
        60,(2*pi*eye(4)/eye(4)).','Dim',4);
    kpoints_name = ["0","2\pi"];

    %EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR_new,"klist_frac",klist_frac);
    Hckt.bandplot((EIGENCAR),OmegaCut,SpectrumX,...
        'title',['4D TI',' n = ',num2str(meshs)],...
        'klist_l',klist_l,'kpoints_l',kpoints_l,'kpoints_name',kpoints_name);
end
%%
end
return;
%%% OpenBoundary %%%%
%% slab
mkdir('slab');
meshs = [20,36,3,3];
% copyfile('KPOINTS_slab','slab/KPOINTS');
% copyfile('POSCAR','slab/POSCAR');
cd('slab');
% Gen_sp
Hckt.Genlib('magnitude',magnitude);
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[1,0,0,0],'analysis',Analysis,'ExciteNodes',[1;2;3;4]);
if strcmp(Analysis,'ac')
        DATAname = [char(Basename),'_AC','.ac0'];
        SPfile = strcat(Basename,'_AC.sp');
        EvalStringDATA = (['!hspiceAC ',DATAname]);
        SPLIS = strcat(Basename,'_AC.lis');
    else
        DATAname = [char(Basename),'.tr0'];
        SPfile = strcat(Basename,'.sp');
        EvalStringDATA = (['!hspiceTR ',DATAname]);
        SPLIS = strcat(Basename,'.lis');
    end
if isunix && ~ismac
% run; run it on linux system
EvalStringRun = (['!hspice -mt 20 -i ',SPfile,' |tee log']);
eval(EvalStringRun);
eval(EvalStringDATA);
else

end
%
if isunix && ~ismac
simulation_result = Hckt.read_hspice_ac(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result,'analysis',Analysis);
%%
% grid
SelectL = VectorList(:,1)==meshs(1) & VectorList(:,end)==1;
% SelectL = VectorList(:,1)==1 & VectorList(:,end)==1;
VectorList2 = VectorList(SelectL,[2,3,4,5]);
ObservationsMat2 = ObservationsMat(SelectL,:);
%
DOSCAR = Hckt.CollectVstruct(VectorList2,ObservationsMat2,OmegaCut,SpectrumX);
[klist_cart,klist_frac,klist_l,kpoints_l,~] = ...
    vasplib.kpathgen(...
    [ ...
     0 1/3 -1/3; ...
     1 1/3 -1/3; ... 
     ],...
    60,(2*pi*eye(3)/eye(3)).','Dim',3);
kpoints_name = ["0","2\pi"];

EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR,"klist_frac",klist_frac,"method","nearest");
Hckt.bandplot(normalize(EIGENCAR),OmegaCut,SpectrumX,...
    'title',['4D TI slab',' n = ',num2str(meshs)],...
    'klist_l',klist_l,'kpoints_l',kpoints_l,'kpoints_name',kpoints_name);
savefig(gcf,'Fig4.fig');
% delete(gcf);
end
%% fermi arc
% Hckt.waveplot(DOSCAR,VectorList2,SpectrumX,[0 0 0],'Frequency',8e06,'Width',0.1e06,'scale',1);
% savefig(gcf,'Fig5.fig');
% delete(gcf);
% save('surf.mat','-v7.3');
% cd(workdir)


%% disk
% mkdir('disk');
% copyfile('KPOINTS_wire','disk/KPOINTS');
% copyfile('POSCAR','disk/POSCAR');
% cd disk;
% % Gen_sp
% Hckt.Genlib();
% [~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[1,1]);
% if strcmp(Analysis,'ac')
%     DATAname = [char(Basename),'_AC','.ac0'];
%     SPfile = strcat(Basename,'_AC.sp');
%     EvalStringDATA = (['!hspiceAC ',DATAname]);
% else
%     DATAname = [char(Basename),'.tr0'];
%     SPfile = strcat(Basename,'.sp');
%     EvalStringDATA = (['!hspiceTR ',DATAname]);
% end
% %
% if isunix && ~ismac
% % run; run it on linux system
% EvalStringRun = (['!hspice -mt 20 -i ',SPfile,' |tee log']);
% eval(EvalStringRun);
% eval(EvalStringDATA);
% end
% %
% if isunix && ~ismac
% simulation_result = Hckt.read_hspice_ac(DATAname);
% [VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result,'analysis',Analysis);
% %
% %% wave plot
% Hckt.waveplot(ObservationsMat,VectorList,SpectrumX,Hcktpre_n.orbL,'Frequency',8e06,'Width',0.1e06,'scale',1);
% savefig(gcf,'Fig5.fig');
% delete(gcf);
% end
cd(workdir)


