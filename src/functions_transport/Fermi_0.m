function f_0 = Fermi_0(Eigs,Ef,T)
% Fermi_Dirac_Distribution
k_B = 8.6173324e-5; % eV/K
tkb = T * k_B;
f_0 = 1./(exp((Eigs - Ef)./tkb) + 1);
f_0(isnan(f_0)) = 0;
end