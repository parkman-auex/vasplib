function orbital_out  = supercell_orb(orbital_init,fin_dir,vacuum_mode)
switch vacuum_mode
    case 0
        orbital_out = orbital_init;
        Ns = [1 0 0;0 1 0;0 0 1];Ns = Ns.*fin_dir;
        fin_dir_list = [0 0 0];
        %disp(fin_dir_list);
        [Rm,sites,Atom_name,Atom_num]=POSCAR_readin();
        [~,~] = supercell(Ns,Rm,sites,Atom_name,Atom_num,fin_dir_list,'POSCAR_super_fin');
        for i = 1:3
            Nslab = fin_dir(i);
            if  Nslab == 0
                Nslab = 1;
            end
            count = 0;
            WAN_NUM = size(orbital_out,1);
            fin_orb = zeros(WAN_NUM*Nslab,3);
            for inum = 1:Nslab     % go over all cells in finite direction
                for j = 1:WAN_NUM  % go over all orbitals in one cell
                    count =count +1;
                    % make a copy of j-th orbital
                    orb_tmp=orbital_out(j,:);
                    % change coordinate along finite direction ; fractional
                    orb_tmp(i)= (orb_tmp(i) + double(inum-1)) / Nslab;
                    % add to the list
                    fin_orb(count,:)=orb_tmp;
                    % do the onsite energies at the same time
                end
            end
            orbital_out = fin_orb;
        end
        
    case 1
        orbital_out = orbital_init;
        for i = 1:3
            Nslab = fin_dir(i);
            if  Nslab == 0
                Nslab = 1;
            end
            count = 0;
            WAN_NUM = size(orbital_out,1);
            fin_orb = zeros(WAN_NUM*Nslab,3);
            for inum = 1:Nslab % go over all cells in finite direction
                for j = 1:WAN_NUM  % go over all orbitals in one cell
                    count =count +1;
                    % make a copy of j-th orbital
                    orb_tmp=orbital_out(j,:);
                    % change coordinate along finite direction ; fractional
                    orb_tmp(i)= (orb_tmp(i) + double(inum-1))/Nslab;
                    % add to the list
                    fin_orb(count,:)=orb_tmp;
                    % do the onsite energies at the same time
                end
            end
            orbital_out = fin_orb;
        end
        
        Ns = [1 0 0;0 1 0;0 0 1];Ns = Ns.*fin_dir;
        fin_dir_list = double(fin_dir>1);
        %disp(fin_dir_list);
        [Rm,sites,Atom_name,Atom_num]=POSCAR_readin();
        % gen POSCAR
        [~,~] = supercell(Ns,Rm,sites,Atom_name,Atom_num,fin_dir_list,'POSCAR_super_fin');
        % rebuild fin_orb
        Rm = Ns*Rm;
        Rmlength1 = norm (Rm(1,:));
        Rmlength2 = norm (Rm(2,:));
        Rmlength3 = norm (Rm(3,:));
        Rm_s_fin_add = [10*Rm(1,:)*fin_dir_list(1)/Rmlength1;...
            10*Rm(2,:)*fin_dir_list(2)/Rmlength2;...
            10*Rm(3,:)*fin_dir_list(3)/Rmlength3];
        Rm_s_fin = Rm + Rm_s_fin_add ;
        Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
        Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
        [nfinorb,~ ]= size(fin_orb);
        for  i = 1:nfinorb
            Rr_orb = fin_orb(i,:)*Rm;
            Rr_s_fin = Rr_orb + Rr_s_fin_add;
            Rc_s_fin = Rr_s_fin / Rm_s_fin;
            fin_orb(i,:) = Rc_s_fin ;
        end
        orbital_out = fin_orb;
end
end