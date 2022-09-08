%% RUN one batch
cd ..
WORKSPACE=pwd;

cd tools ;
POSCAR_dir =  '/Users/parkman/Documents/MATLAB/vasplib/example/example15-1-HTRCI/HT_test/test_one_poscar/C_data/';
C_database_poscar =  strcat(WORKSPACE,'/C_poscar_database/');
C_database_symmetry = strcat(WORKSPACE,'/C_symmetry_database/');
C_database_operation = strcat(WORKSPACE,'/C_operation_database/');
workshop=strcat(WORKSPACE,'/workshop');


phonopy_run = '/Users/parkman/Documents/TOOLs/miniconda2/envs/my_pymatgen/bin/phonopy';
POSCAR_name_struct = dir(POSCAR_dir+"/POSCAR.*");
firstload  = 1;
%

% load old database
if firstload ~= 1
    load('MA_Database.mat');
    count  = length(MA_Database)+1;
    
else % first
    cd(workshop);
    delete('*');
    MA_Database = struct('ID',[],'MA',[]);
    POSCAR_name = POSCAR_name_struct(1).name;
    copyfile( strcat(POSCAR_dir,'/',POSCAR_name), POSCAR_name);
    [POSCAR_file,POSCAR_syminfor_file,ID] = SymInfor_gen(POSCAR_name,phonopy_run);
    MA_Database(1).ID = ID;
    tmp_MA = MaterialAccom(POSCAR_file,POSCAR_syminfor_file,'MaterilProject');
    MA_Database(1).MA = tmp_MA ;
    new_name = tmp_MA.SG+"-"+tmp_MA.SYM_num+"-"+tmp_MA.Atom_num_total+"-"+ID;
    copyfile(POSCAR_file,strcat(C_database_poscar,'/',POSCAR_file));
    copyfile(POSCAR_syminfor_file,strcat(C_database_operation,'/',POSCAR_syminfor_file));
    if tmp_MA.isInversion == 1 
        copyfile(POSCAR_file,strcat(C_database_symmetry,'/',new_name ));
    end
    count  = 2 ;
end
% run POSCAR_name_struct
POSCAR_name_struct_num = length(POSCAR_name_struct );

for i = 1: POSCAR_name_struct_num
    cd(workshop);
    delete('*');
    POSCAR_name = POSCAR_name_struct(i).name;
    copyfile( strcat(POSCAR_dir,'/',POSCAR_name), POSCAR_name);
    [POSCAR_file,POSCAR_syminfor_file,ID] = SymInfor_gen(POSCAR_name,phonopy_run);
    tmp_MA = MaterialAccom(POSCAR_file,POSCAR_syminfor_file,'MaterilProject');
    % MA_Database_num
    MA_Database_num = length(MA_Database);
    %
    %     if MAinMA_Database(tmp_MA,MA_Database)
    %         MA_Database(i).ID = ID;
    %         MA_Database(i).MA = tmp_MA;
    %     end
    add_label = 1;
    for  j =1:MA_Database_num
        if tmp_MA == MA_Database(j).MA
            add_label = 0;
            break;
        end
    end
    % if add label
    if add_label == 1
        MA_Database(count).ID = ID;
        MA_Database(count).MA = tmp_MA;
        new_name = tmp_MA.SG+"-"+tmp_MA.SYM_num+"-"+tmp_MA.Atom_num_total+"-"+ID;
        copyfile(POSCAR_file,strcat(C_database_poscar,'/',POSCAR_file));
        copyfile(POSCAR_syminfor_file,strcat(C_database_operation,'/',POSCAR_syminfor_file));
        if tmp_MA.isInversion == 1 
            copyfile(POSCAR_file,strcat(C_database_symmetry,'/',new_name ));
        end
        count = count+1;
    else
        fprintf('check repeat POSCAR ID: %d\n',ID);
    end
end
MA_Database_num = length(MA_Database);

cd(WORKSPACE);
cd('tools');
if firstload ~= 1
    copyfile('MA_Database.mat','MA_Database.mat.bk')
end
% save('MA_Database.mat','MA_Database');
save('MA_Database.mat','MA_Database');



