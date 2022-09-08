%% RUN_gen_POSCAR;
load('Database_2Dmatpedia.mat');
%% 
Ncrys = length(Database_2Dmatpedia);
target_dir = 'POSCAR_collection';
mkdir(target_dir);

for i = 1:Ncrys
    temp_pymatgen_struct = Database_2Dmatpedia(i).structure;
    ID_string = strsplit(Database_2Dmatpedia(i).material_id,'-');
    ID = str2double(ID_string{2});
    POSCAR_filename = target_dir +"/"+"POSCAR."+string(ID);
    fprintf('(%4d/%-4d)th %s:%s -->%s\n',i,Ncrys,...
        Database_2Dmatpedia(i).material_id,...
        Database_2Dmatpedia(i).chemsys,...
        POSCAR_filename);
    pymatgen2sites(temp_pymatgen_struct,POSCAR_filename );
end