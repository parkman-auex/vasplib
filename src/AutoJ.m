function AutoJ_table = AutoJ(targetdir)
    if nargin < 1
        targetdir = '.';
    end
    workdir = pwd;
    cd(targetdir);
    [Rm,sites,Atom_name,~]=POSCAR_read();
    [~,MAGMOM_string] = system("grep MAGMOM "+ 'INCAR');
    MAGMOM_string = strsplit(MAGMOM_string,'=');
    MAGMOM_str = MAGMOM_string{2};
    MAGMOM_str = strtrim(MAGMOM_str);
    [ncol,Mag_ion_num,sites_mag,sites] = MAGMOM_read(MAGMOM_str,sites);

    search_range = [4, 4, 4];
    Accuracy = 3 ;
    [Atom_store_smart,nn_store_smart,Rnn,Rnn_map]=nn_smart(Rm,sites_mag,search_range,Accuracy,10);

    level_cut  = 4;

    % J_pluslist  = [];
    % J_minuslist = [];
    % J_totlist   = [];
    J_list = zeros(3,level_cut);
    if ncol == 1
        for i = 1:Mag_ion_num
            for j = 1:Mag_ion_num
                bondsign = sites_mag(i).sign * sites_mag(j).sign;
                nn  = nn_store_smart(i,j).nn;
                for k = 1:length(nn)
                    nn_level = nn(k).nn_level;
                    if nn_level <= level_cut
                        switch bondsign 
                            case 1
                                J_list(1,nn_level) =  J_list(1,nn_level) +1 ;
                            case -1
                                J_list(2,nn_level) =  J_list(2,nn_level) - 1 ;
                        end
                    end
                end
            end
        end
        J_list(3,:) = J_list(1,:) + J_list(2,:);
        J_list = J_list/2;    
    elseif ncol == 0 %for more explict classification
        for i = 1:Mag_ion_num
            for j = 1:Mag_ion_num
                bondsign = sites_mag(i).sign * sites_mag(j).sign;
                nn  = nn_store_smart(i,j).nn;
                for k = 1:length(nn)
                    nn_level = nn(k).nn_level;
                    if nn_level <= level_cut
                        switch bondsign
                            case 1
                                J_list(1,nn_level) =  J_list(1,nn_level) +1 ;
                            case -1
                                J_list(2,nn_level) =  J_list(2,nn_level) - 1 ;
                        end
                    end
                end
            end
        end
        J_list(3,:) = J_list(1,:) + J_list(2,:);
        J_list = J_list/2;
    else
        
    end
    [MAG_final_list,U_eV,ENERGY_TOTEN_eV,SYSTEM,NELECT,LSORBIT]= OUTCAR_read('OUTCAR','MAG');
    if ENERGY_TOTEN_eV < 0
        STATUS = "CONVERGE";
    else
        STATUS = "NOT CONVERGE";
    end
    
    sites_final = sites;
    sites_mag_final = sites_mag;
    count = 0;
    NION_sites = length(sites);
    if NION_sites == size(MAG_final_list,1)
        for i = 1:NION_sites
            sites(i).mag_amp = MAG_final_list(i,:);
            if sites(i).mag_ion
                count = count+1;
                sites_mag_final(count).mag_amp = MAG_final_list(i,:);
            end
        end
        cif_gen('POSCAR_init.mcif',Rm,Atom_name,sites,1,'P1_mcif');
        cif_gen('POSCAR_final.mcif',Rm,Atom_name,sites_final,1,'P1_mcif');
        cif_gen('MAG_show_init.mcif',Rm,Atom_name,sites_mag,1,'P1_mcif');
        cif_gen('MAG_show_final.mcif',Rm,Atom_name,sites_mag_final,1,'P1_mcif');
        save('AutoJ.mat');  
    else
        STATUS = "NOT CONVERGE";
    end

    J1 = J_list(3,1);J2 = J_list(3,2);J3 = J_list(3,3);J4 = J_list(3,4);
    NIONS = length(sites);LNONCOLLINEAR = ncol;NMAG_ION = Mag_ion_num;
    ENERGY_REDUCED_eV = ENERGY_TOTEN_eV /  NMAG_ION;

    AutoJ_table = table(SYSTEM,U_eV,ENERGY_TOTEN_eV,NMAG_ION,ENERGY_REDUCED_eV,NIONS,NELECT,LNONCOLLINEAR ,LSORBIT,J1,J2,J3,J4,STATUS);
    cd(workdir);
    
end
