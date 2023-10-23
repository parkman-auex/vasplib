function [ncol,Mag_ion_num,sites_mag,sites] = MAGMOM_read(MAGMOM_str,sites)

    ion_num = length(sites);
    MAGMOM_str = strtrim(MAGMOM_str);
    TEMP_str_list = strsplit(MAGMOM_str);
    Add_num = 0 ;
    Mag_ion_num = 0;
    Mag_init_list = [];
    Mag_seq_list = [];
    count = 0;
    for i =1:length(TEMP_str_list)
        if strcontain(TEMP_str_list{i},"*")
            Temp_str = strsplit(TEMP_str_list{i},'*');
            Add_num = str2double(Temp_str{1})+ Add_num -1;
            count = count + Add_num+1;
        else
            Mag_ion_num = Mag_ion_num +1;
            Mag_init_list =  [Mag_init_list;str2double(TEMP_str_list{i})];
            count = count + 1;
            Mag_seq_list  =   [Mag_seq_list;count];
        end
    end
    tot_num = length(TEMP_str_list) + Add_num ;
    if tot_num/3 == ion_num
        Mag_ion_num = Mag_ion_num /3;
        ncol = 1;
    else
        ncol = 0;
    end
    if ncol 
        if tot_num/3 ~= ion_num
            error('POSCAR and INCAR not match!');
        end
    else
        if tot_num ~= ion_num
            error('POSCAR and INCAR not match!');
        end
    end
    % if
    if ncol == 1
        Mag_init_list = reshape(Mag_init_list,3,Mag_ion_num).';
        Mag_seq_list  =   reshape(Mag_seq_list,3,Mag_ion_num).';
        Mag_seq_list =  Mag_seq_list(:,3)/3;
        count = 0;
        for i = 1:ion_num
            if ismember( i, Mag_seq_list)
                count = count +1;
                sites(i).mag_ion = 1;
                sites(i).mag_amp = Mag_init_list(count,:);
                tmp_sign = sum(sign(sites(count).mag_amp));
                if issame(sites(i).mag_amp,[0,0,0])
                    sites(i).sign = tmp_sign;
                elseif tmp_sign == 0
                    sites(i).sign =sign(sites(i).mag_amp(1))*sign(sites(i).mag_amp(2))^2+...
                        sign(sites(i).mag_amp(2))*sign(sites(i).mag_amp(3))^2+...
                        sign(sites(i).mag_amp(3))*sign(sites(i).mag_amp(1))^2;
                else
                     sites(i).sign = tmp_sign;
                end
    
                sites_mag(count) = sites(i);
                sites_mag(count).seq = count;

                fprintf('(%3d/%-3d) th mag ion: %4s (%5.1f,%5.1f,%5.1f), which is (%3d/%-3d) ion in POSCAR.\n',...
                    count,Mag_ion_num,...
                    sites(i).name,...
                    Mag_init_list(count,1),Mag_init_list(count,2),Mag_init_list(count,3),...
                    i,ion_num);
            else
                sites(i).mag_ion = 0;
            end
        end
    else
        Mag_init_list = reshape(Mag_init_list,1,Mag_ion_num).';
        Mag_seq_list  =   reshape(Mag_seq_list,1,Mag_ion_num).';
        Mag_seq_list =  Mag_seq_list(:,1)/1;
        count = 0;
        for i = 1:ion_num
            if ismember( i, Mag_seq_list)
                count = count +1;
                sites(i).mag_ion = 1;
                sites(i).mag_amp = Mag_init_list(count,:);
                tmp_sign = sum(sign(sites(i).mag_amp));
                sites(i).sign = tmp_sign;
    
                sites_mag(count) = sites(i);
                sites_mag(count).seq = count;

                fprintf('(%3d/%-3d) th mag ion: %4s (%5.1f), which is (%3d/%-3d) ion in POSCAR.\n',...
                    count,Mag_ion_num,...
                    sites(i).name,...
                    Mag_init_list(count,1),...
                    i,ion_num);
            else
                sites(i).mag_ion = 0;
            end
        end
    end
end