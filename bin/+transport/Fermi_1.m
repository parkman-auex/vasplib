function df_dE = Fermi_1(Eigs,Ef,T)
% first derivative of Fermi_Dirac_Distribution 
k_B = 8.6173324e-5; % eV/K
tkb = T * k_B;
%%
ecore = exp((Eigs - Ef)./tkb);
df_dE = -ecore./(tkb .* (ecore + 1).^2);
df_dE(isnan(df_dE)) = 0;
end