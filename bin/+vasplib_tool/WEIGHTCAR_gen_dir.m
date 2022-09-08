function [WEIGHTCAR_struct_cell,Name_list,width] = WEIGHTCAR_gen_dir(mode)
    import vasplib_tool.*;
    directory =dir();
    count  = 0;
    switch mode
        case 'DOS'
            Name_list = [""];
            for i = 1:length(directory)
                filename =  directory(i).name;
                if contains(filename ,'PDOS_') & contains(filename ,'.dat')
                                        count = count+1;
                    Name_list = [Name_list ;strrep(strrep(filename,'.dat',''),'PDOS_','') ];
                    [Pdensity(:,:,count),~,~] = WEIGHTCAR_gen(filename,-1,'vaspkit-DOS-silence');

                end
            end
            Name_list(1,:) = [];
            WEIGHTCAR_struct_cell = Pdensity;
            [~,width] = size(Pdensity(:,:,1));
        case 'PBAND'
            Name_list = [""];
            for i = 1:length(directory)
                filename =  directory(i).name;
                if contains(filename ,'PBAND_') & contains(filename ,'.dat')
                    count = count+1;
                    Name_list = [Name_list ;strrep(strrep(filename,'.dat',''),'PBAND_','') ];
                    [WEIGHTCAR_3d_mat{count},~,~] = WEIGHTCAR_gen(filename,-1,'vaspkit-band-silence');

                end
            end
            Name_list(1,:) = [];
            width= size(WEIGHTCAR_3d_mat{1},3);

            for i = 1:count
                WEIGHTCAR_struct_cell{i} = temp_3d_mat2struct(WEIGHTCAR_3d_mat{i},Name_list(i,:),width);
            end
    end
            
end

function WEIGHTCAR_struct = temp_3d_mat2struct(WEIGHTCAR_3d,Name,width)
    import vasplib_tool.*
    % f or not
    if width >10
        num2orbital_name =  orbital_maprule_gen(1);
    else
        num2orbital_name = orbital_maprule_gen(0);
    end
    for i = 1:width
        WEIGHTCAR_struct(i).WEIGHTCAR = WEIGHTCAR_3d(:,:,i);
        WEIGHTCAR_struct(i).displayname = Name + "-"+num2orbital_name(i);
    end
end