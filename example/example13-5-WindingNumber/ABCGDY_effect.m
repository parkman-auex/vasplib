%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%% 
syms Delta m v lambda u1 u2;
syms k_x k_y k_z d real ;
H_2D = HK(4,2);
%% 
% % the coeffs is not allowed use x y z

H_2D = H_2D ...
    + Term( -Delta + m*k_x^2+m*k_y^2 ,    sigma_z*tau_0)...
    + Term( v*k_x,                   sigma_x*tau_x)...
    + Term( v*k_y,                   sigma_x*tau_y)...
    + Term( lambda*(k_x^2-k_y^2),    sigma_z*tau_x)...
    + Term( 2*lambda*(k_x * k_y),    sigma_z*tau_y)...
    ;
H_3D = H_2D ...
    + Trig( -u1*(cos(d*k_z) )         ,    sigma_z*tau_0)...
    + Trig( -u2*(sin(d*k_z))          ,    sigma_y*tau_z)...
    ;
H_3D = H_3D < 'POSCAR';
H_3D = H_3D < 'KPOINTS';
%% 
% % para

b = 1.424e-10;
h_bar =6.582119514e-16; %eV⋅s
v = 7e5 * h_bar /b  ; % m/s

%% 
d =3.3;
% 
%v =0;
m = 0.21;
lambda = 0.00; % eV
u1 =0.233;
u2 =0.26;
Delta =0.213;
tolerance = 0.01;
%kmesh = [10,10,0.825525678548721];
H_3Dn = H_3D.Subsall();
%
EIGENCAR = H_3Dn.EIGENCAR_gen();
bandplot(EIGENCAR,[-1,1]);

%[k_list_cart,klist_f,~]=H_3Dn.findnodes(kmesh,2,tolerance);
%%
[nodes_s, nodes_r] = findnodes(H_3Dn,'nk',[8,8,8],'original_point',[-0.5,-0.5,-0.5],'Num_Occupied',2);
%%
H_3Dn.BZplot('color','none')
scatter3(nodes_r(:,1),nodes_r(:,2),nodes_r(:,3),5,'filled');
%
GammaOper = double(-sigma_x*tau_z);
% kloop= H_3Dn.kloop_gen([0 0.0593 0.070707;1 0 1;0.00593 0 101],'kplane');
% kloopz = linspace(0,0,101).';
% kloop = [zeros(101,2),kloopz ];
% kloop = kloop * H_3Dn.Gk;
kpoint_r = [0,0.0327017,0.824835];
[klist_loop_r,klist_loop] = vasplib.kloop1D(kpoint_r,[1,0,0],0.01,'inputCar',true,'enforceCar',true);
plot3(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3),'Color','k')
BF = H_3Dn.BP_1D(klist_loop);
fprintf('Berryphase along loop:%7.5f.\n',BF);
%
[klist_loop_r,~] = vasplib.kloop1D(kpoint_r,[1,0,0],0.02,'inputCar',true,'enforceCar',true,'nk',501);
plotg(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3));
[WindingNumber,~] = H_3Dn.WindingNumber(GammaOper,klist_loop_r);
fprintf('Upper ring: Winding Number along loop:%7.5f.\n',WindingNumber);
axis off;
% axis equal
kpoint_r = [0,0.0327017,-0.824835];
[klist_loop_r,~] = vasplib.kloop1D(kpoint_r,[1,0,0],0.02,'inputCar',true,'enforceCar',true,'nk',501);
plotg(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3));
[WindingNumber,~] = H_3Dn.WindingNumber(GammaOper,klist_loop_r);
fprintf('Upper ring: Winding Number along loop:%7.5f.\n',WindingNumber);
% POSCAR_read;
%%
[kloop,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(H_3Dn.Rm,'KPOINTS_winding');
kpoints_r = [0.0374837921933913,0.0216412775131593,0.134625947812235];
A = TopoCharge(H_3Dn,GammaOper,kloop);