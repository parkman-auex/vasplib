%% Declare the primitive of QTI
TITLE = 'QTI';
Dim = 2;
meshs = [10,10];
magnitude = 'u';
workdir = pwd;
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
TBmodelForHckt_QTI;
%% TRANSLATE
HCKT_test = Hcktpre_n.HR2Hckt('title',TITLE,'dim',Dim,'autominus',true,'magnitude',magnitude);
%%
if 1 == 1
%% Gen_sp
Hckt.Genlib('magnitude',magnitude)
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs);
DATAname = [char(Basename),'.tr0'];
end
if 1 == 1
%% run; run it on linux system
EvalStringRun = (['!hspice ',strcat(Basename,'.sp'),' |tee log']);
eval(EvalStringRun);
EvalStringDATA = (['!hspiceTR ',DATAname]);
eval(EvalStringDATA);
end
%% 
if 1 ==1
%% data clean
DATAname = [Basename,'.tr0'];
simulation_result = Hckt.read_hspice_tr(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result);
%% grid
SelectL = VectorList(:,end) == 1;
ObservationsMat2 = ObservationsMat(SelectL,:);
[DOSCAR_3D,~,~,OmgL] = Hckt.CollectVstruct2D(VectorList,ObservationsMat2,OmegaCut,SpectrumX);
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR_3D);
Hckt.bandplot(EIGENCAR,OmegaCut,OmgL);
%%
end

%%% OpenBoundary %%%%
%% slab
cd slab
% Gen_sp
Hckt.Genlib('magnitude',magnitude);
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[0,1]);
DATAname = [char(Basename),'.tr0'];
%
if 1 == 1
% run; run it on linux system
EvalStringRun = (['!hspice ',strcat(Basename,'.sp'),' |tee log']);
eval(EvalStringRun);
EvalStringDATA = (['!hspiceTR ',DATAname]);
eval(EvalStringDATA);
end
%
if 1==1
DATAname = [Basename,'.tr0'];
simulation_result = Hckt.read_hspice_tr(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result);
%
%% grid
SelectL = VectorList(:,end)== 4 & VectorList(:,end-1)== 10 ;
% VectorList(:,2)== 30;
ObservationsMat2 = ObservationsMat(SelectL,:);
%[~,ObservationsMat2] = HollowKnight.generalcontractrow2(VectorList(SelectL,2),ObservationsMat2);
OmgL = SpectrumX;
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
if 1 == 1
% run; run it on linux system
EvalStringRun = (['!hspice ',strcat(Basename,'.sp'),' |tee log']);
eval(EvalStringRun);
EvalStringDATA = (['!hspiceTR ',DATAname]);
eval(EvalStringDATA);
end
%
if 1==1
DATAname = [Basename,'.tr0'];
simulation_result = Hckt.read_hspice_tr(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result);
%
%% grid
SelectNode = [10,10];
SelectL = VectorList(:,end)== 1 & VectorList(:,1) == SelectNode(1) & VectorList(:,2) == SelectNode(2);
% VectorList(:,2)== 30;
ObservationsMat2 = ObservationsMat(SelectL,:);
Hckt.Frequencyplot(abs(ObservationsMat2),SpectrumX,OmegaCut);
end
cd(workdir)


