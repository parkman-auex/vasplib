%% **************** 1D **************** 
syms t k_x real;
H_1D = Htrig(1,Trig(t*cos(k_x),eye(1)),'Dim',1);
[H_1D.klist_cart,H_1D.klist_frac,klist_l,kpoints_l,~] = ...
    vasplib.kpathgen([-0.5;0;0;0.5],60,H_1D.Gk,'Dim',1);
t = -1;
H_1D_n = H_1D.Subsall;
EIGENCAR = H_1D_n.EIGENCAR_gen();
vasplib_plot.bandplot(EIGENCAR,[-1,1],klist_l,kpoints_l,["-\pi","0","\pi"],'Color','r','title','Dim = 1 (Htrig)');
%%
H_1D.orbL = [0];
H_1D_hr = H_1D.Htrig2HR();
% H_1D_hr_list = H_1D_hr.rewrite();
H_1D_hr.list;
H_1D_hr_n = H_1D_hr.Subsall;
EIGENCAR = H_1D_hr_n.EIGENCAR_gen();
vasplib_plot.bandplot(EIGENCAR,[-1,1],klist_l,kpoints_l,["-\pi","0","\pi"],'Color','r','title','Dim = 1 (HR)');

%% **************** 2D **************** 
syms t lambda k_x k_y real;
useful_matrics("sigma",'mode','sym');

H_2D_hr = HR(4,'Dim',2);
H_2D_hr.Rm = [1 0;-1/2 sqrt(3)/2];
H_2D_hr.orbL = [1/3 2/3;2/3 1/3;1/3 2/3;2/3 1/3;];

V_NN = [0 0;1 0;0 -1];
V_NNN = [1 0;0 1;-1 -1];
% t
H_2D_hr = H_2D_hr.set_hop(...
     t,2,1,V_NN,'sym');
H_2D_hr = H_2D_hr.set_hop(...
     t,4,3,V_NN,'sym');
% SOC
H_2D_hr = H_2D_hr.set_hop(...
     kron(sigma_z, -1i*lambda),1,1,V_NNN,'sym');
H_2D_hr = H_2D_hr.set_hop(...
     kron(sigma_z,  1i*lambda),2,2,V_NNN,'sym');
% auto hermi
H_2D_hr = H_2D_hr.autohermi();

% kpath
[H_2D_hr.klist_cart,H_2D_hr.klist_frac,klist_l,kpoints_l,~] = ...
    vasplib.kpathgen(...
    [0 0;0.5 0.0;...
    0.5 0.0;1/3 1/3;...
    1/3 1/3;0 0],...
    60,H_2D_hr.Gk,'Dim',2);
% para and plot
t = -1;
lambda = 0.06;
H_KaneMele_n = H_2D_hr.Subsall;
EIGENCAR = H_KaneMele_n.EIGENCAR_gen();
vasplib_plot.bandplot(EIGENCAR,[-3,3],klist_l,kpoints_l,["\Gamma","M","K","\Gamma"],'Color','r','title','Dim = 2 (HR)');
%%
H_2D = H_2D_hr.HR2Htrig();
simplify(sym(H_2D))
%%
H_2D_hk = H_2D_hr.HR2HK([1/3 1/3],"Order",2);
simplify(sym(H_2D_hk))

%% **************** 3D **************** 
syms t lambda k_x k_y k_z real;
useful_matrics("sigma",'mode','sym');

H_3D = HR(4,'Dim',3);
H_3D.Rm = [sqrt(3)/2 -1/2 0;sqrt(3)/2 1/2 0;0 0 1];
H_3D.orbL = [0 0 0;1/3 1/3 0;0 0 0;1/3 1/3 0;];

V_NN = [0 0 0;1 0 0;0 1 0];
V_NNN = [-1 0 0;0 1 0;1 -1 0];
% t
H_3D = H_3D.set_hop(...
     t,2,1,V_NN,'sym');
H_3D = H_3D.set_hop(...
     t,4,3,V_NN,'sym');
% SOC
H_3D = H_3D.set_hop(...
     kron(sigma_z, -1i*lambda),1,1,V_NNN,'sym');
H_3D = H_3D.set_hop(...
     kron(sigma_z,  1i*lambda),2,2,V_NNN,'sym');
% auto hermi
H_3D = H_3D.autohermi();

% kpath
[H_3D.klist_cart,H_3D.klist_frac,klist_l,kpoints_l,~] = ...
    vasplib.kpathgen(...
    [0 0 0;0.0 0.5 0;...
    0.0 0.5 0;1/3 2/3 0;...
    1/3 2/3 0;0 0 0],...
    60,H_3D.Gk,'Dim',3);
% para and plot
t = -1;
lambda = 0.06;
H_KaneMele_n = H_3D.Subsall;
EIGENCAR = H_KaneMele_n.EIGENCAR_gen();
vasplib_plot.bandplot(EIGENCAR,[-3,3],klist_l,kpoints_l,["\Gamma","M","K","\Gamma"],'Color','r','title','Dim = 3 (HR)');

%% **************** 4D **************** 