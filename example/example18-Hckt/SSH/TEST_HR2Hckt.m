%% Declare the primitive of Octagraphene
TITLE = 'SSH';
Dim = 1;
meshs = [10];
workdir = pwd;
%magnitude = 'p';
magnitude = 'p';
Analysis = 'tran';
switch magnitude
    case 'p'
        OmegaCut = [0 2e7];
    case 'n'
        OmegaCut = [0 6e5];
    case 'u'
        OmegaCut = [0 2e4];
end

%magnitude = 'u';
% E_onsite = Subckt('X1 n+_1 TOGND E_A VarC0=Cc\n');
% E_v = Subsckt('X12 n+_1 n+_2 TOGND C\n');
%%  run TB part
TBmodelForHckt_SSH;
%% TRANSLATE
HCKT_test = Hcktpre_n.HR2Hckt('title',TITLE,'dim',Dim,'autominus',true,'magnitude',magnitude);
%%
if 1 == 1
%% Gen_sp
Hckt.Genlib('magnitude',magnitude);
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs);
DATAname = [char(Basename),'.tr0'];
end
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
% grid
SelectL = VectorList(:,end) == 1;
ObservationsMat2 = ObservationsMat(SelectL,:);
SelectOmegaL = find(SpectrumX >= OmegaCut(1) & SpectrumX <= OmegaCut(2));
OmgL = SpectrumX(SelectOmegaL);
DOSCAR_2D = Hckt.CollectVstruct1D(ObservationsMat2);
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR_2D);
Hckt.bandplot(EIGENCAR,OmegaCut,OmgL);
%%
end

return;

%% AC
if 1 == 1
% Gen_sp
Analysis = 'ac';
Hckt.Genlib('magnitude',magnitude);
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'analysis',Analysis);
DATAname = [char(Basename),'.ac0'];
end
if 1 == 2
%% run; run it on linux system
EvalStringRun = (['!hspice ',strcat(Basename,'.sp'),' |tee log']);
eval(EvalStringRun);
EvalStringDATA = (['!hspiceAC ',DATAname]);
eval(EvalStringDATA);
end
%% 
if 1 ==2
%% data clean
Basename = 'SSH_10';
DATAname = [Basename,'.ac0'];
simulation_result = Hckt.read_hspice_ac([Basename,'.ac0']);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result,'analysis','ac');
% grid
SelectL = VectorList(:,2)==1;
ObservationsMat2 = ObservationsMat(SelectL,:);
%% 
DOSCAR = Hckt.CollectVstruct1D((ObservationsMat2));
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR);
Hckt.bandplot(EIGENCAR,OmegaCut,SpectrumX,'title',['SSH AC n = ',num2str(meshs)]);

end
%%
%%% OpenBoundary %%%%
%% slab
mkdir('slab');
cd slab
% Gen_sp
Hckt.Genlib();
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[1]);
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
end
cd ..;
%%% OpenBoundary AC %%%%
%% slab
mkdir('slab');
cd slab
% Gen_sp
Hckt.Genlib();
[~,Basename] = HCKT_test.hspice_gen('mesh',meshs,'fin_dir',[1],'analysis','ac');
DATAname = [char(Basename),'.ac0'];
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
DATAname = [Basename,'.ac0'];
simulation_result = Hckt.read_hspice_tr(DATAname);
[VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result);
%
end
