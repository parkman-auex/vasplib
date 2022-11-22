function useful_constants()

%% constants
h_eV_s = 4.1357e-15; % eV.s
hbar_eV_s = 6.5821e-16; % eV.s
charge_C = 1.6e-19; % Coulomb
muB_eV_T = 5.7884e-5; % eV/Tesla

%% conventions

%%
assignin('base',"h_eV_s",    h_eV_s);
assignin('base',"hbar_eV_s", hbar_eV_s);
assignin('base',"charge_C",  charge_C);
assignin('base',"muB_eV_T",  muB_eV_T);
end