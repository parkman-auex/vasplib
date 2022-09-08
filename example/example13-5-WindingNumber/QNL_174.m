%% Quadratic and cubic nodal lines stabilized by crystalline symmetry
% https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.121106 (5)
% 10.1103/PhysRevLett.110.240404 (2013)
% 10.1103/PhysRevB.78.195125 (2008)
% Heff(k) = a * k_-^2+ H.c.,
% \mathcal{N}=\frac{1}{4 \pi i} \oint_{C} \operatorname{Tr} \sigma_{z} \mathcal{H}_{\mathrm{eff}}^{-1}(\boldsymbol{q}) \nabla_{\boldsymbol{q}} \mathcal{H}_{\mathrm{eff}}(\boldsymbol{q}) \cdot d \boldsymbol{q}
%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
sigma_plus = double(sigma_x) + 1i*double(sigma_y);
sigma_minus = double(sigma_x) - 1i*double(sigma_y);
%% 
syms a u1 u2 real;
syms k_x k_y k_z d real ;
k_minus = (k_x - 1i*k_y) ;
k_plus = (k_x + 1i*k_y) ;
H_174 = HK(2,2);
%%
H_174 = H_174 ...
    + Term( a*k_minus^2 ,    sigma_plus)...
    + Term( a*k_plus^2  ,    sigma_minus)...
    ;
% H_174  =H_174  ...
%     + Trig( u1*(cos(d*k_z) )         ,    sigma_x)...
%     + Trig( u2*(sin(d*k_z) )         ,    sigma_y)...
%     ;
%%
H_174 = H_174  < 'POSCAR';
H_174 = H_174  < 'KPOINTS';
d =3.3;
a = 0.2;
% u1 =0.233;
% u2 =0.002;
u1 =0;
u2= 0;
tolerance = 0.01;
H_174n = H_174.Subsall();
%
EIGENCAR = H_174n.EIGENCAR_gen();
bandplot(EIGENCAR,[-1,1]);
%%
% kzf = 0:0.01:1;
% H_174n  = H_174n < 'KPOINTS_findnode';
% [k_list_cart,klist_f,gaplist]=H_174n.findnodes2(kzf,1,tolerance);
%%
 [nodes_s, nodes_r] = findnodes(H_174n,'nk',[8,8,8],'original_point',[-0.5,-0.5,-0.5],'Num_Occupied',1);
 H_174n.BZplot('color','none')
 scatter3(nodes_r(:,1),nodes_r(:,2),nodes_r(:,3),5,'filled');
% return;
%%
%H_174n.BZplot();
kpoint_r = [0,0,0];
[klist_loop_r,klist_loop] = vasplib.kloop1D(kpoint_r,[0,0,1],0.1,'inputCar',true,'enforceCar',true);
plot3(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3),'Color','k')
BF = H_174n.BP_1D(klist_loop);
fprintf('Berryphase along loop:%7.5f.\n',BF);
%
[klist_loop_r,klist_loop] = vasplib.kloop1D(kpoint_r,[0,0,1],0.15,'inputCar',true,'enforceCar',true);
plotg(klist_loop_r(:,1),klist_loop_r(:,2),klist_loop_r(:,3));
GammaOper = double(sigma_z);
[WindingNumber,WL] = H_174n.WindingNumber(GammaOper,klist_loop_r);
fprintf('Winding Number along loop:%7.5f.\n',WindingNumber);

axis equal