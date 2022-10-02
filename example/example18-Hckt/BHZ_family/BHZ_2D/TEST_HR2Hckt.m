%% Declare the primitive of BHZ_2D
TITLE = 'BHZ_2D';
Dim = 2;
meshs = [20,20];
magnitude = 'p';
workdir = pwd;
Analysis = 'tran';
switch magnitude
    case 'p'
        OmegaCut = [0 2e7];
    case 'n'
        OmegaCut = [0 6e5];
    case 'u'
        OmegaCut = [0 2e4];
    otherwise
        OmegaCut = [0 2e7];
end
%%  run TB part
TBmodelForHckt_BHZ_2D;
%% TRANSLATE
HCKT_test = Hcktpre_n.HR2Hckt('title',TITLE,'dim',Dim,'autominus',true,'magnitude',magnitude);
%%
if 1 == 1
    %% Gen_sp
    Hckt.Genlib('magnitude',magnitude);
    [~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'analysis',Analysis);
    if strcmp(Analysis,'ac')
        DATAname = [char(Basename),'_AC','.ac0'];
        SPfile = strcat(Basename,'_AC.sp');
        EvalStringDATA = (['!hspiceAC ',DATAname]);
    else
        DATAname = [char(Basename),'.tr0'];
        SPfile = strcat(Basename,'.sp');
        EvalStringDATA = (['!hspiceTR ',DATAname]);
    end
end
if isunix && ~ismac()
%% run; run it on linux system
EvalStringRun = (['!hspice -mt 20 -i ',SPfile,' |tee log']);
eval(EvalStringRun);
eval(EvalStringDATA);
end
%% 
if isunix && ~ismac
%% data clean
simulation_result = Hckt.read_hspice_ac(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result,'analysis',Analysis);
%% grid
SelectL = VectorList(:,3)==1;
ObservationsMat2 = ObservationsMat(SelectL,:);
VectorList2 = VectorList(SelectL,:);
%% 
DOSCAR = Hckt.CollectVstruct2D(VectorList2,ObservationsMat2,OmegaCut,SpectrumX);
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR);
Hckt.bandplot(EIGENCAR,OmegaCut,SpectrumX,'title',[TITLE,' n = ',num2str(meshs)]);
end

%%% OpenBoundary %%%%
%% slab
mkdir('slab');
copyfile('KPOINTS_slab','slab/KPOINTS');
copyfile('POSCAR','slab/POSCAR')
% Gen_sp
Hckt.Genlib('magnitude',magnitude);
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[0,1]);
DATAname = [char(Basename),'.tr0'];
%
% Gen_sp
Hckt.Genlib();
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[0,1],'analysis',Analysis);
if strcmp(Analysis,'ac')
    DATAname = [char(Basename),'_AC','.ac0'];
    SPfile = strcat(Basename,'_AC.sp');
    EvalStringDATA = (['!hspiceAC ',DATAname]);
else
    DATAname = [char(Basename),'.tr0'];
    SPfile = strcat(Basename,'.sp');
    EvalStringDATA = (['!hspiceTR ',DATAname]);
end
if isunix && ~ismac
% run; run it on linux system
EvalStringRun = (['!hspice -mt 20 -i ',SPfile,' |tee log']);
eval(EvalStringRun);
eval(EvalStringDATA);
end
%
if isunix && ~ismac
simulation_result = Hckt.read_hspice_ac(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result,'analysis',Analysis);
%
% grid
SelectL = VectorList(:,2)==meshs(1) & VectorList(:,3)==1;
ObservationsMat2 = ObservationsMat(SelectL,:);
%% 
DOSCAR = Hckt.CollectVstruct1D((ObservationsMat2));
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR);
Hckt.bandplot(EIGENCAR,OmegaCut,SpectrumX,'title','slab');

end
cd(workdir)
%% disk
mkdir('disk');
copyfile('KPOINTS_wire','disk/KPOINTS');
copyfile('POSCAR','disk/POSCAR');
% Gen_sp
Hckt.Genlib();
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[1,1]);
if strcmp(Analysis,'ac')
    DATAname = [char(Basename),'_AC','.ac0'];
    SPfile = strcat(Basename,'_AC.sp');
    EvalStringDATA = (['!hspiceAC ',DATAname]);
else
    DATAname = [char(Basename),'.tr0'];
    SPfile = strcat(Basename,'.sp');
    EvalStringDATA = (['!hspiceTR ',DATAname]);
end
%
if isunix && ~ismac
% run; run it on linux system
EvalStringRun = (['!hspice -mt 20 -i ',SPfile,' |tee log']);
eval(EvalStringRun);
eval(EvalStringDATA);
end
%
if isunix && ~ismac
simulation_result = Hckt.read_hspice_ac(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result,'analysis',Analysis);
%
%% wave plot
Hckt.waveplot(ObservationsMat,VectorList,SpectrumL,Hcktpre_n.orbL,'Frequency',7.17e06,'Width',0.5e06,'scale',8);
end
cd(workdir)


