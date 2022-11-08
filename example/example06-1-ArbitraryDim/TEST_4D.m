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

% gamma_0 = sigma_0 * tau_0 ;
% gamma_1 = sigma_x * tau_x ;
% gamma_2 = sigma_x * tau_y ;
% gamma_3 = sigma_x * tau_z ; 
% gamma_4 = sigma_y * tau_0 ;
% gamma_5 = sigma_z * tau_0 ;

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
     0 0 1/3 -1/3; ...
     0 1 1/3 -1/3; ...
     0 0 0 0; ...
     0 0 1 0; ...
     0 0 0 0; ...
     0 0 0 1; ...
     ],...
    60,H_4D.Gk,'Dim',4);
kpoints_name = ["\Gamma","\Gamma_x|\Gamma","\Gamma_y|\Gamma","\Gamma_z|\Gamma","\Gamma_w"];

% −t/2<m<t,C2 =−2
%
%% para
t = 1;
m = 0.0;
epsilon = 9;

H_4D_n = H_4D.Subsall();

% plot

EIGENCAR = H_4D_n.EIGENCAR_gen();
vasplib_plot.bandplot(EIGENCAR,[3,14],klist_l,kpoints_l,kpoints_name,'Color','r','title','Dim = 4 (Htrig)');

%%
H_4D.Rm = eye(4);
H_4D.orbL = zeros(4);
H_4D_hr = H_4D.Htrig2HR();
H_4D_hr_n = H_4D_hr.Subsall();

EIGENCAR = H_4D_hr_n.EIGENCAR_gen();
vasplib_plot.bandplot(EIGENCAR,[3,14],klist_l,kpoints_l,kpoints_name,'Color','r','title','Dim = 4 (HR)');

%% slab
% we use supercell_orb and supercell_hr for the arbitrary Dim
H_4D_hr_slab_n = H_4D_hr_n.supercell_hr(diag([40,1,1,1]),'OBC',[1 0 0 0]);
% H_4D_hr_slab_n = H_4D_hr_n.rewrite.supercell_hr(diag([40,1,1,1])); % also check list mode
%% kpath
[H_4D_hr_slab_n.klist_cart,H_4D_hr_slab_n.klist_frac,klist_l,kpoints_l,~] = ...
    vasplib.kpathgen(...
    [ ...
     0 0 0 0; ...
     1 0 0 0; ... 
     0 0 1/3 -1/3; ...
     0 1 1/3 -1/3; ...
     0 5/12 0 -1/3; ...
     0 5/12 1 -1/3; ...
     0 5/12 1/3 0; ...
     0 5/12 1/3 1; ...
     ],...
    60,H_4D.Gk,'Dim',4);
kpoints_name = ["\Gamma","\Gamma_x|\Gamma","\Gamma_y|\Gamma","\Gamma_z|\Gamma","\Gamma_w"];
% plot

EIGENCAR = H_4D_hr_slab_n.EIGENCAR_gen();
vasplib_plot.bandplot(EIGENCAR,[3,14],klist_l,kpoints_l,kpoints_name,'Color','r','title','Dim = 4 -slab');

