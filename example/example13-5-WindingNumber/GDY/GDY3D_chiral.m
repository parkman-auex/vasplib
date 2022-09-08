%% Parity TB test for GDY3
% * 进行 nn_sk 搜索
% * 这里 search_range 取 【1 1 0】，Accuracy 取小数点3位，maxR值取 5 埃

search_range = [1 1 1];
Accuracy = 3;
r_max_search = 6;
GDY3D = HR.from_POSCAR_SE("POSCAR",...
    'Type','mat',...
    'r_max',r_max_search,...
    'level_cut',30,...
    'search_range',search_range,...
    'Accuracy',1e-3,...
    'chiral',true,...
    'onsite',false,'WAN_NUM',18);
%% 
% * 查看一下体系的未知量（symbolic）
% * GDY3D.symvar_list
%% 参数
Rnn = GDY3D.nn_information;
VppP_1 = 3.34;
VppP_2 = 2.99;
VppP_3 = 2.82;
VppP_4 = 2.71*0.6;

VppP_5 = 0.0;
VppP_6 = 0;
VppP_7 = 0;
VppP_8 = 0;
VppP_9 = 0;
VppP_11 = 0;

VppP_10 = 0.1;
VppP_12 = 0.1;


VppS_1 = 0;
VppS_2 = 0;
VppS_3 = 0;
VppS_4 = 0;
VppS_5 = 0;
VppS_6 = 0;
VppS_7 = 0;
VppS_8 = 0;
VppS_9 = 0;
VppS_10 = 0;

VppS_11 = 0.05*1.4;
VppS_12 = 0.42*1.4;

VppP_0  = 2.7;
VppS_0  = -0.48*1.7;

%Rnn
a0 = 1.42 ;delta0=0.45255; % VppP
d0 = 3.35; 
%delta0=0;
for i = 13:32
    VppP_n = "VppP_"+string(i);
    VppS_n = "VppS_"+string(i);
    tmp_string =  VppP_n+"  = VppP_0*exp(-(Rnn(i)-a0)/delta0) ";
    tmp_string2 =   VppS_n+" = VppS_0*exp(-(Rnn(i)-d0)/delta0)  ";
    eval(tmp_string);
    eval(tmp_string2);
end
%% 能带

GDY3D_n = Subsall(GDY3D);
% bulk band
%%
GDY3D_n = GDY3D_n <'KPOINTS';
EIGENCAR =GDY3D_n.EIGENCAR_gen();
% plot
bandplot(EIGENCAR,[-0.5,0.5]);
%% 寻找nodeline
[nodes_s, nodes_r] = findnodes(GDY3D_n,'nk',[8,8,8],'original_point',[-0.5,-0.5,-0.5],'Num_Occupied',9);
%%
GDY3D_n.BZplot('color','none')
scatter3(nodes_r(:,1),nodes_r(:,2),nodes_r(:,3),5,'filled');
return;
%%
GammaOper = kron([1 0;0 -1],eye(9));
% kloop= H_3Dn.kloop_gen([0 0.0593 0.070707;1 0 1;0.00593 0 101],'kplane');
% kloopz = linspace(0,0,101).';
% kloop = [zeros(101,2),kloopz ];
% kloop = kloop * H_3Dn.Gk;
kpoint_r = [-0.000011789113372  -0.0000000   0.852768102575870];
[klist_loop_r,klist_loop] = vasplib.kloop1D(kpoint_r,[0,1,0],1e-5,'inputCar',true,'enforceCar',true);
plot3(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3),'Color','k')
BF = GDY3D_n.BP_1D(klist_loop);
fprintf('Berryphase along loop:%7.5f.\n',BF);
%
[klist_loop_r,~] = vasplib.kloop1D(kpoint_r,[0,1,0],1e-5,'inputCar',true,'enforceCar',true,'nk',501);
plotg(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3));
[WindingNumber,~] = GDY3D_n.WindingNumber(GammaOper,klist_loop_r);
fprintf('Upper ring: Winding Number along loop:%7.5f.\n',WindingNumber);
axis off;
% axis equal
kpoint_r = [-0.000011789113372  -0.0000000   0.852768102575870];
[klist_loop_r,~] = vasplib.kloop1D(kpoint_r,[0,1,0],1e-5,'inputCar',true,'enforceCar',true,'nk',501);
plotg(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3));
[WindingNumber,~] = GDY3D_n.WindingNumber(GammaOper,klist_loop_r);
fprintf('Lower ring: Winding Number along loop:%7.5f.\n',WindingNumber);
% kzf = 0.0  :0.01:1;
% GDY3D_n = GDY3D_n < 'KPOINTS_findnode';
% [k_list_cart,klist_f,gaplist,fig]=GDY3D_n.findnodes2(kzf);
%% 拓扑数
% HRz_gen

% GammaOper = kron([1 0;0 -1],eye(9));
% % kloop= GDY3D_n.kloop_gen([0.013,  0.013, 0.435;1 0 0;0.0065 0 101],'kplane');
% kloopz = linspace(0,1,101).';
% kloop = [zeros(101,2),kloopz ];
% kloop = kloop * GDY3D_n.Gk;
% kpoints_r = [0.00881061059071106	0.0152604251888630,0.884525848855331];
% GDY3D_n.TopoCharge(GammaOper,kloop,kpoints_r);