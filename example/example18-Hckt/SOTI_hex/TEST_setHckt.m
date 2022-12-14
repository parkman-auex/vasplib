%% SOTIHEX component
%% declare a Hckt
HCKT_SOTIHEX = Hckt('title','SOTIHEX','Nports',6,'vectorAll',[0,0]);
%% set Lib
HCKT_SOTIHEX.Lib = [HCKT_SOTIHEX.Lib ;];
%% set options
HCKT_SOTIHEX.Options = [".ic v(n_01_01_4) = 1";...
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
HCKT_SOTIHEX = HCKT_SOTIHEX.set_hop([-1,0],CW,4,1);
HCKT_SOTIHEX = HCKT_SOTIHEX.set_hop([0,-1],CW,5,2);
HCKT_SOTIHEX = HCKT_SOTIHEX.set_hop([-1,-1],CW,6,3);
%% Gen_sp
meshs = [10,10];
% for i = [10,20,50,100]
%     meshs=i;
Hckt.Genlib();
meshname = ['_',num2str(meshs(1)),'_',num2str(meshs(2))];
   [~,Basename] = HCKT_SOTIHEX.hspice_gen(...
        "probenode",'Allnode',...
        'mesh',meshs,'mode','vectorized');
    DATAname = [char(Basename),'.tr0'];
% end
%% 
if isunix && ~ismac()
    %% run; run it on linux system
    EvalStringRun = (['!hspice -mt 20 -i ',strcat(Basename,'.sp'),' |tee log']);
    eval(EvalStringRun);
    EvalStringDATA = (['!hspiceTR ',DATAname]);
    eval(EvalStringDATA);
    meshs = [10,10];
    meshname = ['_',num2str(meshs(1)),'_',num2str(meshs(2))];
    simulation_result = Hckt.read_hspice_tr(DATAname);
    [VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result);
%ObservationsMat2 = ObservationsMat;
%% grid
L_0 = 1e-6;
Cfactor = 100*1E-12;% 100 pF
OmegaFactor = (Cfactor*L_0*100)^(-1/2);
OMEGACUTpre = [0,2];
OmegaCut = OMEGACUTpre*OmegaFactor;
SelectL = VectorList(:,3)== 2;
ObservationsMat2 = ObservationsMat(SelectL,:);
[DOSCAR_3D,klist1,klist2,OmgL] = Hckt.CollectVstruct2D(VectorList,ObservationsMat2,OmegaCut,SpectrumX);
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR_3D);
Hckt.bandplot(EIGENCAR,OmegaCut,OmgL);
end
return;

%%% AC
%%
   [~,Basename] = HCKT_SOTIHEX.hspice_gen(...
        "probenode",'Allnode',...
        'analysis','ac','mesh',meshs);
    DATAname = [char(Basename),'_AC.ac0'];
 if isunix && ~ismac()
    %% run; run it on linux system
    EvalStringRun = (['!hspice -mt 20 -i ',strcat(Basename,'_AC.sp'),' |tee log']);
    eval(EvalStringRun);
    EvalStringDATA = (['!hspiceAC ',DATAname]);
    eval(EvalStringDATA);
    meshs = [10,10];
    meshname = ['_',num2str(meshs(1)),'_',num2str(meshs(2))];
    simulation_result = Hckt.read_hspice_ac(DATAname);
    [VectorList,ObservationsMat,~,SpectrumX,TimeL] = Hckt.extractObservations(simulation_result,'analysis','ac');
%ObservationsMat2 = ObservationsMat;
%% grid
L_0 = 1e-6;
Cfactor = 100*1E-12;% 100 pF
OmegaFactor = (Cfactor*L_0*100)^(-1/2);
OMEGACUTpre = [0,2];
OmegaCut = OMEGACUTpre*OmegaFactor;
SelectL = VectorList(:,3)== 2;
ObservationsMat2 = ObservationsMat(SelectL,:);
[DOSCAR_3D,klist1,klist2,OmgL] = Hckt.CollectVstruct2D(VectorList,ObservationsMat2,OmegaCut,SpectrumX);
EIGENCAR = Hckt.ProjectDOSCAR(DOSCAR_3D);
Hckt.bandplot(EIGENCAR,OmegaCut,OmgL);
end   