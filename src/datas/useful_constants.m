function useful_constants()

%% constants
h_eV_s = 4.1357e-15; % eV.s
hbar_eV_s = 6.5821e-16; % eV.s
charge_C = 1.6e-19; % Coulomb
muB_eV_T = 5.7884e-5; % eV/Tesla
kB_eV_K = 8.61733e-5; % eV/Kelvin

assignin('base',"h_eV_s",    h_eV_s);
assignin('base',"hbar_eV_s", hbar_eV_s);
assignin('base',"charge_C",  charge_C);
assignin('base',"muB_eV_T",  muB_eV_T);
assignin('base',"kB_eV_K",   kB_eV_K);
%% conventions

%% matlab_default_colors
blue_dark = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
purple = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
blue_light = [0.3010 0.7450 0.9330];
red = [0.6350 0.0780 0.1840];

assignin('base',"blue_dark", blue_dark);
assignin('base',"orange", orange);
assignin('base',"yellow", yellow);
assignin('base',"purple", purple);
assignin('base',"green", green);
assignin('base',"blue_light", blue_light);
assignin('base',"red", red);
%%


end