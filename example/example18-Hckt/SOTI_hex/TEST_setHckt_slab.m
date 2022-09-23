%% SOTIHEX component
%% declare a Hckt
HCKT_SOTIHEX = Hckt('title','SOTIHEX','Nports',6,'vectorAll',[0,0]);
%% set Lib
HCKT_SOTIHEX.Lib = [HCKT_SOTIHEX.Lib ;];
%% set options
HCKT_SOTIHEX.Options = ["*.ic v(n_01_01_1) = 1";...
    ".tran 1ns 10us";...
    ".option post=2 probe";...
    ".option parhier=global";...
    ".option accurate=1";...
    "*.option INTERP";...
    ".PARAM InitV =0V VarC0 =100p C_hopping = 100p C_w = 200p C_v = 100p"];
%% set_homecell
% homecell
HomeCell = Subckt([ ...
    'X1 1 TOGND E_A\n',...
    'Cv1 1 3 C_v\n',...
    'X2 2 TOGND E_A\n' ...
    'Cv2 3 2 C_v\n',...
    'X3 3 TOGND E_A\n',...
    'Cv3 2 4 C_v\n' ...
    'X4 4 TOGND E_A\n',...
    'Cv4 4 6 C_v\n',...
    'X5 5 TOGND E_A\n', ...
    'Cv5 6 5 C_v\n', ...
    'X6 6 TOGND E_A\n', ...
    'Cv6 5 1 C_v\n', ...
    ],...
    'Pri',...
    '1 2 3 4 5 6 TOGND',...
    'VarC0 =100p VarL0 = 1u InitV = 0V C_v = 100p' ...
    );
% Subcktobj,PortInL,PortOutL,DescriptionL
HCKT_SOTIHEX = HCKT_SOTIHEX.set_home(HomeCell,[1,2,3],[4,5,6]);
%% set_hop
CW = Subckt('CBA IN N0 C_w');
% vector,Subcktobj,PortInL,PortOutL,DescriptionL
HCKT_SOTIHEX = HCKT_SOTIHEX.set_hop([1,0],CW,1,4);
HCKT_SOTIHEX = HCKT_SOTIHEX.set_hop([0,1],CW,2,5);
HCKT_SOTIHEX = HCKT_SOTIHEX.set_hop([1,1],CW,3,6);
%% Gen_sp
meshs = [10,10];
% for i = [10,20,50,100]
%     meshs=i;
basename = 'SOTIhex';
meshname = ['_',num2str(meshs(1)),'_',num2str(meshs(2))];
NAME = [basename,meshname];
Hckt.Genlib();
HCKT_SOTIHEX.hspice_gen([NAME,'.sp'],...
    "probenode",'Allnode',...
    'mesh',meshs,'mode','general','fin_dir',[0,1],'node',2);
% end

if isunix && ~ismac()
    %%
    %eval(['!hspice ',NAME,'.sp >log']);
    %eval(['!hspiceTR ',NAME,'.tr0']);
    %ObservationsMat2 = ObservationsMat;
    EvalStringRun = (['!hspice -mt 20 -i ',strcat(NAME,'.sp'),' |tee log']);
    eval(EvalStringRun);
    DATAname = [NAME,'.tr0'];
    EvalStringDATA = (['!hspiceTR ',DATAname]);
    eval(EvalStringDATA);
    meshs = [10,10];
    basename = 'SOTIhex';
    meshname = ['_',num2str(meshs(1)),'_',num2str(meshs(2))];
    NAME = [basename,meshname]; 
    simulation_result = Hckt.read_hspice_tr([NAME,'.tr0']);
    [VectorList,ObservationsMat,~,OmgL,TimeL] = Hckt.extractObservations(simulation_result);

%% grid
L_0 = 1e-6;
Cfactor = 100*1E-12;% 100 pF
OmegaFactor = (Cfactor*L_0*100)^(-1/2);
OMEGACUTpre = [0,2];
OmegaCut = OMEGACUTpre*OmegaFactor;
 SelectL = VectorList(:,3)== 2 &VectorList(:,2)== meshs(2);
% VectorList(:,2)== 30;
ObservationsMat2 = ObservationsMat(SelectL,:);
DOSCAR = Hckt.CollectVstruct1D(ObservationsMat2);
%[DOSCAR_3D,klist1,klist2,OmgL] = Hckt.CollectVstruct2D(VectorList,ObservationsMat2,OmegaCut,SpectrumX);
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR,'KPOINTS','KPOINTS_slab');
Hckt.bandplot(EIGENCAR,OmegaCut,OmgL,'KPOINTS','KPOINTS_slab');
end

%return;

%%% AC
%%

HCKT_SOTIHEX.hspice_gen([NAME,'_AC.sp'],...
    "probenode",'Allnode',...
    'mesh',meshs,'mode','general','analysis','ac','fin_dir',[0,1],'node',2);
% end

if isunix && ~ismac()
    %%
    %eval(['!hspice ',NAME,'.sp >log']);
    %eval(['!hspiceTR ',NAME,'.tr0']);
    EvalStringRun = (['!hspice -mt 20 -i ',strcat(NAME,'_AC.sp'),' |tee log']);
    eval(EvalStringRun);
    DATAname = [NAME,'_AC.ac0'];
    EvalStringDATA = (['!hspiceAC ',DATAname]);
    eval(EvalStringDATA);
    meshs = [10,10];
    basename = 'SOTIhex';
    meshname = ['_',num2str(meshs(1)),'_',num2str(meshs(2))];
    NAME = [basename,meshname];
    simulation_result = Hckt.read_hspice_ac([NAME,'_AC.ac0']);
    [VectorList,ObservationsMat,~,OmgL,TimeL] = Hckt.extractObservations(simulation_result,'analysis','ac');
    %ObservationsMat2 = ObservationsMat;
    
    %% grid
    L_0 = 1e-6;
    Cfactor = 100*1E-12;% 100 pF
    OmegaFactor = (Cfactor*L_0*100)^(-1/2);
    OMEGACUTpre = [0,2];
    OmegaCut = OMEGACUTpre*OmegaFactor;
    SelectL = VectorList(:,3)== 2 &VectorList(:,2)== meshs(2) ;
    % VectorList(:,2)== 30;
    ObservationsMat2 = ObservationsMat(SelectL,:);
  
    DOSCAR = Hckt.CollectVstruct1D(ObservationsMat2);
    %[DOSCAR_3D,klist1,klist2,OmgL] = Hckt.CollectVstruct2D(VectorList,ObservationsMat2,OmegaCut,SpectrumX);
    EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR,'KPOINTS','KPOINTS_slab');
    Hckt.bandplot(EIGENCAR,OmegaCut,OmgL,'KPOINTS','KPOINTS_slab');
end