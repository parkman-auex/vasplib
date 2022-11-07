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
H_K = simplify(sym(H_2D_hk));
syms k_p k_m real;
simplify(subs(H_K,[k_x k_y],[(k_p+k_m)/2 (k_p-k_m)/2i]))

%% **************** 3D **************** 
% syms t lambda k_x k_y k_z real;
% useful_matrics("sigma",'mode','sym');
% H_3D = HR(4,'Dim',3);
% H_3D.Rm = [sqrt(3)/2 -1/2 0;sqrt(3)/2 1/2 0;0 0 1];
% H_3D.orbL = [0 0 0;1/3 1/3 0;0 0 0;1/3 1/3 0;];
% 
% V_NN = [0 0 0;1 0 0;0 1 0];
% V_NNN = [-1 0 0;0 1 0;1 -1 0];
% % t
% H_3D = H_3D.set_hop(...
%      t,2,1,V_NN,'sym');
% H_3D = H_3D.set_hop(...
%      t,4,3,V_NN,'sym');
% % SOC
% H_3D = H_3D.set_hop(...
%      kron(sigma_z, -1i*lambda),1,1,V_NNN,'sym');
% H_3D = H_3D.set_hop(...
%      kron(sigma_z,  1i*lambda),2,2,V_NNN,'sym');
% % auto hermi
% H_3D = H_3D.autohermi();
% 
% % kpath
% [H_3D.klist_cart,H_3D.klist_frac,klist_l,kpoints_l,~] = ...
%     vasplib.kpathgen(...
%     [0 0 0;0.0 0.5 0;...
%     0.0 0.5 0;1/3 2/3 0;...
%     1/3 2/3 0;0 0 0],...
%     60,H_3D.Gk,'Dim',3);
% % para and plot
% t = -1;
% lambda = 0.06;
% H_KaneMele_n = H_3D.Subsall;
% EIGENCAR = H_KaneMele_n.EIGENCAR_gen();
% vasplib_plot.bandplot(EIGENCAR,[-3,3],klist_l,kpoints_l,["\Gamma","M","K","\Gamma"],'Color','r','title','Dim = 3 (HR)');

%% **************** 4D **************** 
%  10.1093/nsr/nwaa065
% \mathcal{H}(\boldsymbol{k})=\sum_{a=0}^5 f_a(\boldsymbol{k}) \gamma_a
% γ_1, 2, 3 = τ_1, 2, 3 ⊗ρ1, γ_4 = τ_0⊗ρ2, and γ_5 = τ_0⊗ρ3
useful_matrics(["sigma","tau"]);

gamma_0 = tau_0 * sigma_0;
gamma_1 = tau_x * sigma_x;
gamma_2 = tau_y * sigma_x;
gamma_3 = tau_z * sigma_x;
gamma_4 = tau_0 * sigma_y;
gamma_5 = tau_0 * sigma_z;

% f0(k) = ε − t cos(k2 + k3),
% f1(k) = −t(1 + cos k1 + cos k2),
% f2(k) = t(sink1 + sink2), 
% f3(k) = −t(1 + cosk3 + cos k4), 
% f4(k) = t(sin k3 + sin k4), 
% f5(k) = m − t cos(k2 + k3),

syms epsilon t m k_x k_y k_z k_w real;

f0k = epsilon - t*cos(k_y + k_z)  ;
f1k = -t*(1 + cos(k_x) + cos(k_y));
f2k =  t*(sin(k_x) + sin(k_y))    ;   
f3k = -t*(1 + cos(k_z) + cos(k_w)); 
f4k =  t*(sin(k_z) + sin(k_w))    ;
f5k =       m - t*cos(k_y + k_z)  ;


H_4D = Htrig(4,'Dim',4);

H_4D = H_4D ...
   +  Trig(f0k,gamma_0)...
   +  Trig(f1k,gamma_1)...
   +  Trig(f2k,gamma_2)...
   +  Trig(f3k,gamma_3)...
   +  Trig(f4k,gamma_4)...
   +  Trig(f5k,gamma_5)...
;

% kpath
[H_4D.klist_cart,H_4D.klist_frac,klist_l,kpoints_l,~] = ...
    vasplib.kpathgen(...
    [ ...
     0 0 0 0; ...
     1 0 0 0; ... 
     0 0 0 0; ...
     0 1 0 0; ...
     0 0 0 0; ...
     0 0 1 0; ...
     0 0 0 0; ...
     0 0 0 1; ...
     ],...
    60,H_4D.Gk,'Dim',4);
kpoints_name = ["\Gamma","\Gamma_x|\Gamma","\Gamma_y|\Gamma","\Gamma_z|\Gamma","\Gamma_w"];

% −t/2<m<t,C2 =−2
%
% para
t = 1;
m = 0;
epsilon = 0;

H_4D_n = H_4D.Subsall();

% plot

EIGENCAR = H_4D_n.EIGENCAR_gen();
vasplib_plot.bandplot(EIGENCAR,[-6,5],klist_l,kpoints_l,kpoints_name,'Color','r','title','Dim = 4 (Htrig)');
