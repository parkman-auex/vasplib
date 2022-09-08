function [POSCAR_file,POSCAR_syminfor_file,ID] = SymInfor_gen(POSCAR_name,phonopy_run)
%% [POSCAR_file,POSCAR_syminfor_file,ID] = SymInfor_gen(POSCAR_name,phonopy_run)
if nargin <2
    phonopy_run = '/Users/parkman/Documents/TOOLs/miniconda2/envs/my_pymatgen/bin/phonopy';
end

eval("!"+phonopy_run+' --tolerance 0.01 --symmetry -c '+POSCAR_name+...
    '>'+'SymInfor_'+POSCAR_name);
try
    copyfile('PPOSCAR','POSCAR');
catch
    disp('wrong POSCAR');
    ID = -1;
    POSCAR_file = '';
    POSCAR_syminfor_file = '';
    return;
end
eval("!"+phonopy_run+' --tolerance 0.01 --symmetry -c '+'POSCAR'+...
    '>'+'SymInfor_'+POSCAR_name);
delete('BPOSCAR');
delete('POSCAR');
try
    movefile('PPOSCAR',POSCAR_name);
catch
    disp('wrong POSCAR');
    
    return;
end

POSCAR_file = POSCAR_name;
POSCAR_syminfor_file = "SymInfor_"+POSCAR_name ;
ID          = str2double(strrep(POSCAR_file,'POSCAR.',''));
end