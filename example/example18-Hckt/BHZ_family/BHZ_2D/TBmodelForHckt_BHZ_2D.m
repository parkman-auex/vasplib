%% armchair 2
%% BHZ model
% Useful matrix

% clear;
useful_matrics;

syms M0 M1 M2 real;
syms A B real;
syms k_x k_y k_z real;
%
M       = M0-M2*(k_x^2+k_y^2)-M1*(k_z^2);
k_plus  = k_x + 1i* k_y;
k_minus = k_x - 1i* k_y;
%
BHZ = HK(4,2);
%
BHZ = BHZ ...
    +Term(-A*k_x ,sigma_y*tau_x )...
    +Term(-A*k_y ,sigma_0*tau_y )...
    +Term(M     ,sigma_0*tau_z )...
    +Term(B*(k_x^2-k_y^2),sigma_x*tau_x )...
    ;
% Bulk band
% para
BHZ = BHZ <'POSCAR';
BHZ_TB= BHZ.kp2TB();
%%

HcktTB = BHZ_TB.rewrite;
%%

%% Para

%%
% unit eV 
M0     =  1   ;
%M0     = -0.25   ;
% unit eV Ang
A      =  1      ;
B      =  0    ;
% unit eV Ang^2
M1     =  0      ;
M2     =  0.5  ;
% unit    Ang
a      =  1        ;
b      =  1       ;
c      =  1        ;
% Zhao prb
C_0 = 3;
%% prepare
Hcktpre = HcktTB.HRforHckt();
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


%%
FigTest = Figs(1,2);


% drawnow;
bandplot(HcktTB.Subsall().EIGENCAR_gen('convention','I'),...
    [-5,5],'Color','b','ax',FigTest.axes(1),'ylabel','Energy');
%
Hcktpre_n = Hcktpre.rewind.Subsall();
%
EIGENCAR = Hcktpre_n.EIGENCAR_gen();
F0CAR = (EIGENCAR*L_0*Cfactor).^(-1/2)./(2*pi);
OMEGACUTpre = [0,2];
OMEGACUT = OMEGACUTpre*OmegaFactor;
%
bandplot(F0CAR,OMEGACUTpre*OmegaFactor,'ylabel','Frequency (Hz)','ax',FigTest.axes(2));
%%
FigTest = Figs(1,2);
repeatnum   = 20;
fin_dir     =  2;
[EIGENCAR_slab,klist_l,kpoints_l,kpoints_name] = Hcktpre_n.slab(repeatnum,fin_dir,'KPOINTS_slab');
% plot
F0CAR_slab = (EIGENCAR_slab*L_0*Cfactor).^(-1/2)./(2*pi);
bandplot(F0CAR_slab,OMEGACUTpre*OmegaFactor,klist_l,kpoints_l,kpoints_name,'ylabel','Energy','title','slab','ax',FigTest.axes(1));
%%
Nsuper = 20;
Nslab = [Nsuper Nsuper 1];
n_disk = Hcktpre_n.Hnanowire_gen(Nslab);
[EIGENCAR_disk,WAVECAR_disk,WEIGHTCAR] = n_disk.EIGENCAR_gen('WEIGHTCAR',true,'klist',[0 0 0]);
F0CAR_disk = (EIGENCAR_disk*L_0*Cfactor).^(-1/2)./(2*pi);
scatter(FigTest.axes(2),1:length(EIGENCAR_disk),F0CAR_disk,100*ones(length(EIGENCAR_disk),1),(WEIGHTCAR).^3,'filled');
colormap(ColorMap.Matplotlib('coolwarm'))
Nocc = length(EIGENCAR_disk)/2;
xlim([Nocc-Nsuper+1 Nocc+Nsuper]);