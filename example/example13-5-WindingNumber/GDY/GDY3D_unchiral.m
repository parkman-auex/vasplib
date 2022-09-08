HRobj = HR.from_wannier90('wannier90_hr.dat');
HRobj.input_Rm('POSCAR');

%% 能带
GDY3D_n = HRobj;
% bulk band
%%
GDY3D_n = GDY3D_n <'KPOINTS';
EIGENCAR =GDY3D_n.EIGENCAR_gen();
% plot
bandplot(EIGENCAR,[-0.5,0.5]);
%% 寻找nodeline
% [nodes_s, nodes_r] = findnodes(GDY3D_n,'nk',[8,8,8],'original_point',[-0.5,-0.5,-0.5],'Num_Occupied',9);
NodesDATA = importdata('Nodes.dat');
nodes_r = NodesDATA.data;
GDY3D_n.BZplot('color','none')
scatter3(nodes_r(:,1),nodes_r(:,2),nodes_r(:,3),5,'filled');
return;
%%
kpoint_r = [-0.643676,0.351277,-0.159846];
[klist_loop,klist_loop_r] = vasplib.kloop1D(kpoint_r,[1,1,0],0.01,'enforceCar',true);
plot3(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3),'Color','k');
BF = HRobj.BP_1D(klist_loop);
fprintf('Berryphase along loop:%7.5f.\n',BF);
%%
GammaOper = kron([1 0;0 -1],eye(36));

[klist_loop,klist_loop_r] = vasplib.kloop1D(kpoint_r,[1,1,0],0.01+0.002,'enforceCar',true);
plotg(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3));
[WindingNumber,~] = GDY3D_n.WindingNumber(GammaOper,klist_loop_r);
fprintf('Upper ring: Winding Number along loop:%7.5f.\n',WindingNumber);
axis off;
% axis equal
kpoint_r = [0.643676,-0.351277,0.159846];
[klist_loop,klist_loop_r] = vasplib.kloop1D(kpoint_r,[1,1,0],0.01+0.002,'enforceCar',true);
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