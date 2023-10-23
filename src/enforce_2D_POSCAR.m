%%
function enforce_2D_POSCAR(origin_POSCAR,finial_POSCAR,vacuum_length,minum_bond_length)

    [Rm,sites,Atom_name,Atom_num,~]=POSCAR_readin(origin_POSCAR);
    [~,~,Rnn,~]=nn_smart(Rm,sites,[1 1 1],3);
    minum_bond = Rnn(1);
    if nargin < 4
        minum_bond_length =  minum_bond;
    end
    Rm = Rm*minum_bond_length/minum_bond;
    % check c
    Rf_list = [[sites.rc1]',[sites.rc2]',[sites.rc3]'];

    var_a = var(Rf_list(:,1));
    var_b = var(Rf_list(:,2));
    var_c = var(Rf_list(:,3));
    
    var_a_length = norm(var_a*Rm(1,:));
    var_b_length = norm(var_b*Rm(2,:));
    var_c_length = norm(var_c*Rm(3,:));
    
    Rm_length1 = norm(Rm(1,:));
    Rm_length2 = norm(Rm(2,:));
    Rm_length3 = norm(Rm(3,:));
    
    var_min = min([var_a,var_b,var_c]);
    
    test_label = 0;
    
    switch var_min
        case var_a
            c_dir = 1;
            if Rm_length1 >= vacuum_length
                fprintf('exchange a with c')
                test_label = 1;
                C_dir = 1;
            end
        case var_b
            c_dir = 2;
            if Rm_length2 >= vacuum_length
                fprintf('exchange b with c')
                test_label = 1;
                C_dir = 2;
            end
        case var_c
            c_dir = 3;
            if Rm_length3 >= vacuum_length
                fprintf('no need to change')
                test_label = 1;
                C_dir = 3;
            end
        otherwise
            test_label = 2;
    end
    if test_label == 0 
        [~,c_dir] = max([Rm_length1,Rm_length2,Rm_length3]);
        fprintf('use max length in %d dir',c_dir);
    end
    if test_label == 1|test_label == 0 
        switch c_dir
            case 1
                if Rm_length1 == max([Rm_length1,Rm_length2,Rm_length3])
                    Rm(1,:) = Rm(1,:)*vacuum_length/Rm_length1;
                    C_dir = 1;
                else
                    fprintf('Program cant do anything, take your own!');
                    return;
                end
            case 2
                if Rm_length2 == max([Rm_length1,Rm_length2,Rm_length3])
                    Rm(2,:) = Rm(2,:)*vacuum_length/Rm_length2;
                    C_dir = 2;
                else
                    fprintf('Program cant do anything, take your own!');
                                        return;
                end
            case 3                
                if Rm_length3 == max([Rm_length1,Rm_length2,Rm_length3])
                    Rm(3,:) = Rm(3,:)*vacuum_length/Rm_length3;
                    C_dir = 3;
                else
                    fprintf('Program cant do anything, take your own!');
                                    return;
                end
        end
        
    end
    %     P_matr
    % exchange
    switch C_dir
        case 3
            POSCAR_gen(Rm,sites,Atom_name,Atom_num,finial_POSCAR);
        case 2
            P_matrix = [1 0 0;0 0 1;0 1 0];
            Rm = P_matrix * Rm * P_matrix.';
            Rf_list = Rf_list * P_matrix ;
            for i =1:length(sites)
                sites(i).rc1 =  Rf_list(i,1);
                sites(i).rc2 =  Rf_list(i,2);
                sites(i).rc3 =  Rf_list(i,3);
            end
            POSCAR_gen(Rm,sites,Atom_name,Atom_num,finial_POSCAR);
        case 1
            P_matrix = [0 0 1;0 1 0;1 0 0];
            Rm = P_matrix * Rm * P_matrix.';
            Rf_list = Rf_list * P_matrix ;
            for i =1:length(sites)
                sites(i).rc1 =  Rf_list(i,1);
                sites(i).rc2 =  Rf_list(i,2);
                sites(i).rc3 =  Rf_list(i,3);
            end
            POSCAR_gen(Rm,sites,Atom_name,Atom_num,finial_POSCAR); 
    end
    
end