%% Declare the primitive of Octagraphene
TITLE = 'QTI';
Dim = 2;
meshs = [2,2];
OmegaCut = [0 2e7];
workdir = pwd;
magnitude = 'u';
switch magnitude
    case 'p'
        OmegaCut = [0 2e7];
    case 'n'
        OmegaCut = [0 6e5];
    case 'u'
        OmegaCut = [0 2e4];
end
% E_onsite = Subckt('X1 n+_1 TOGND E_A VarC0=Cc\n');
% E_v = Subsckt('X12 n+_1 n+_2 TOGND C\n');
%%  run TB part
TBmodelForHckt_QTI2;
%% TRANSLATE
HCKT_test = Hcktpre_n.HR2Hckt('title',TITLE,'dim',Dim,'autominus',true,'magnitude',magnitude);
%%
if 1 == 1
%% Gen_sp
Hckt.Genlib('magnitude',magnitude);
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs);
DATAname = [char(Basename),'.tr0'];
end
error();
if 1 == 2
%% run; run it on linux system
EvalStringRun = (['!hspice ',strcat(Basename,'.sp'),' |tee log']);
eval(EvalStringRun);
EvalStringDATA = (['!hspiceTR ',DATAname]);
eval(EvalStringDATA);
end
%% 
if 1 ==2
%% data clean
DATAname = [Basename,'.tr0'];
simulation_result = Hckt.read_hspice_tr(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result);
%% grid
SelectL = VectorList(:,end) == 4;
ObservationsMat2 = ObservationsMat(SelectL,:);
[DOSCAR_4D,~,~,~,OmgL] = Hckt.CollectVstruct3D(VectorList,ObservationsMat2,OmegaCut,SpectrumX);
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR_4D);
Hckt.bandplot(EIGENCAR,OmegaCut,OmgL);
%%
end

return;
%%% OpenBoundary %%%%
%% slab
cd slab
% Gen_sp
Hckt.Genlib();
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[0,1]);
DATAname = [char(Basename),'.tr0'];
%
if 1 == 2
% run; run it on linux system
EvalStringRun = (['!hspice ',strcat(Basename,'.sp'),' |tee log']);
eval(EvalStringRun);
EvalStringDATA = (['!hspiceTR ',DATAname]);
eval(EvalStringDATA);
end
%
if 1==2
DATAname = [Basename,'.tr0'];
simulation_result = Hckt.read_hspice_tr(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result);
%
%% grid
SelectL = VectorList(:,end)== 4;
% VectorList(:,2)== 30;
ObservationsMat2 = ObservationsMat(SelectL,:);
ObservationsMat2 = HollowKnight.generalcontractrow(VectorList(SelectL,2),ObservationsMat2);
DOSCAR = Hckt.CollectVstruct1D(ObservationsMat2);
%[DOSCAR_3D,klist1,klist2,OmgL] = Hckt.CollectVstruct2D(VectorList,ObservationsMat2,OmegaCut,SpectrumX);
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR);
Hckt.bandplot(EIGENCAR,OmegaCut,OmgL);
end
cd(workdir)
%% SOTI
cd disk
% Gen_sp
Hckt.Genlib();
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[1,1]);
DATAname = [char(Basename),'.tr0'];
%
if 1 == 2
% run; run it on linux system
EvalStringRun = (['!hspice ',strcat(Basename,'.sp'),' |tee log']);
eval(EvalStringRun);
EvalStringDATA = (['!hspiceTR ',DATAname]);
eval(EvalStringDATA);
end
%
if 1==2
DATAname = [Basename,'.tr0'];
simulation_result = Hckt.read_hspice_tr(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result);
%
%% grid
SelectL = VectorList(:,end)== 4;
% VectorList(:,2)== 30;
ObservationsMat2 = ObservationsMat(SelectL,:);
end
cd(workdir)
