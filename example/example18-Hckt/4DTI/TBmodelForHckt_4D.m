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
%
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
H_4D.Rm = eye(4);
H_4D.orbL = zeros(4);
HcktTB = H_4D.Htrig2HR();
%% kpath
[HcktTB.klist_cart,HcktTB.klist_frac,klist_l,kpoints_l,~] = ...
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
    60,HcktTB.Gk,'Dim',4);
kpoints_name = ["\Gamma","\Gamma_x|\Gamma","\Gamma_y|\Gamma","\Gamma_z|\Gamma","\Gamma_w"];
%% prepare
Hcktpre = HcktTB.HRforHckt("coefficient",1);
%% Add Unit
L_0 = 1e-6;
MU = 100;
try
    switch magnitude
        case 'p'
            Cfactor = MU*1E-12;% 100 pF
        case 'n'
            Cfactor = MU*1E-09;% 100 nF
        case 'u'
            Cfactor = MU*1E-06;% 100 uF
        case 'm'
            Cfactor = MU*1E-03;% 100 mF
    end
catch
    Cfactor = MU*1E-06;% 100 uF
end
OmegaFactor = (Cfactor*L_0*100)^(-1/2);
t = -1;
m = 0.0;
C_0 = 0;
epsilon = 9;
%%
FigTest = Figs(1,2);

Hcktpre_n = Hcktpre.Subsall();
%
EIGENCAR = Hcktpre_n.EIGENCAR_gen();
F0CAR = (EIGENCAR*L_0*Cfactor).^(-1/2)./(2*pi);
OMEGACUTpre = [0.4,0.8];
OMEGACUT = OMEGACUTpre*OmegaFactor;
%
bandplot(F0CAR,OMEGACUTpre*OmegaFactor,klist_l,kpoints_l,kpoints_name,'ylabel','Frequency (Hz)','ax',FigTest.axes(1));

%%
H_4D_hr_slab_n = Hcktpre_n.supercell_hr(diag([30,1,1,1]),'OBC',[1 0 0 0]);
EIGENCAR_slab = H_4D_hr_slab_n.EIGENCAR_gen();
% plot
OMEGACUTpre = [0.4,0.7];
F0CAR_slab = (EIGENCAR_slab*L_0*Cfactor).^(-1/2)./(2*pi);
bandplot(F0CAR_slab,OMEGACUTpre*OmegaFactor,klist_l,kpoints_l,kpoints_name,'ylabel','Energy','title','slab');
savefig(gcf,'Fig2.fig');
delete(gcf);