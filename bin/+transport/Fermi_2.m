function d2f_dE2 = Fermi_2(Eigs,Ef,T)
% second derivative of Fermi_Dirac_Distribution 
k_B = 8.6173324e-5; % eV/K
tkb = T * k_B;
%%
ecore = exp(-(Ef - Eigs)./tkb);
d2f_dE2 = (2.*exp(-2.*(Ef - Eigs)./tkb))./(tkb^2.*(ecore + 1).^3) - ecore./(tkb^2.*(ecore + 1).^2);
d2f_dE2(isnan(d2f_dE2)) = 0;
end