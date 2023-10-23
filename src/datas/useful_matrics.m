%% useful_tools
function useful_matrics(InputStrlist,options)
arguments
    InputStrlist = ["sigma","tau","s","pauli","Gamma","Gell-Mann"];
    options.mode {mustBeMember(options.mode,{'double','sym','pauli_matric'})} = 'pauli_matric';
end
F = eval(['@',options.mode]);

sigma_0 = F(pauli_matric(0));
sigma_x = F(pauli_matric(1));
sigma_y = F(pauli_matric(2));
sigma_z = F(pauli_matric(3));

Gamma_0 = F(gamma_matric(0));
Gamma_1 = F(gamma_matric(1));
Gamma_2 = F(gamma_matric(2));
Gamma_3 = F(gamma_matric(3));
Gamma_4 = F(gamma_matric(4));
Gamma_5 = F(gamma_matric(5));
Gamma_12 = F(gamma_matric(1,2));
Gamma_13 = F(gamma_matric(1,3));
Gamma_14 = F(gamma_matric(1,4));
Gamma_15 = F(gamma_matric(1,5));
Gamma_23 = F(gamma_matric(2,3));
Gamma_24 = F(gamma_matric(2,4));
Gamma_25 = F(gamma_matric(2,5));
Gamma_34 = F(gamma_matric(3,4));
Gamma_35 = F(gamma_matric(3,5));
Gamma_45 = F(gamma_matric(4,5));

Lambda_0 = F(Gell_Mann(0));
Lambda_1 = F(Gell_Mann(1));
Lambda_2 = F(Gell_Mann(2));
Lambda_3 = F(Gell_Mann(3));
Lambda_4 = F(Gell_Mann(4));
Lambda_5 = F(Gell_Mann(5));
Lambda_6 = F(Gell_Mann(6));
Lambda_7 = F(Gell_Mann(7));
Lambda_8 = F(Gell_Mann(8));

if sum(contains(InputStrlist,"sigma"|"pauli"))
    assignin('base',"sigma_0",sigma_0);
    assignin('base',"sigma_x",sigma_x);
    assignin('base',"sigma_y",sigma_y);
    assignin('base',"sigma_z",sigma_z);
end

if sum(contains(InputStrlist,"tau"))
    assignin('base',"tau_0",sigma_0);
    assignin('base',"tau_x",sigma_x);
    assignin('base',"tau_y",sigma_y);
    assignin('base',"tau_z",sigma_z);
end

if sum(contains(InputStrlist,"s"))
    assignin('base',"s_0",sigma_0);
    assignin('base',"s_x",sigma_x);
    assignin('base',"s_y",sigma_y);
    assignin('base',"s_z",sigma_z);
end

if sum(contains(InputStrlist,"Gamma"))
    assignin('base',"Gamma_0",Gamma_0);
    assignin('base',"Gamma_1",Gamma_1);
    assignin('base',"Gamma_2",Gamma_2);
    assignin('base',"Gamma_3",Gamma_3);
    assignin('base',"Gamma_4",Gamma_4);
    assignin('base',"Gamma_5",Gamma_5);

    assignin('base',"Gamma_12",Gamma_12);
    assignin('base',"Gamma_13",Gamma_13);
    assignin('base',"Gamma_14",Gamma_14);
    assignin('base',"Gamma_15",Gamma_15);

    assignin('base',"Gamma_23",Gamma_23);
    assignin('base',"Gamma_24",Gamma_24);
    assignin('base',"Gamma_25",Gamma_25);

    assignin('base',"Gamma_34",Gamma_34);
    assignin('base',"Gamma_35",Gamma_35);

    assignin('base',"Gamma_45",Gamma_45);
end

if sum(contains(InputStrlist,"Gell-Mann"))
    assignin('base',"Lambda_0",Lambda_0);
    assignin('base',"Lambda_1",Lambda_1);
    assignin('base',"Lambda_2",Lambda_2);
    assignin('base',"Lambda_3",Lambda_3);
    assignin('base',"Lambda_4",Lambda_4);
    assignin('base',"Lambda_5",Lambda_5);
    assignin('base',"Lambda_6",Lambda_6);
    assignin('base',"Lambda_7",Lambda_7);
    assignin('base',"Lambda_8",Lambda_8);
end

end