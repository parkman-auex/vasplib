%% make sure you have copy the POSCAR KPOINTS EIGENVAL DOSCAR in this directory
ls;
if isfile('POSCAR')
     fprintf('POSCAR<ok>')% File exists.
else
     fprintf('POSCAR<error>')% File does not exist.
end
if isfile('KPOINTS')
     fprintf('KPOINTS<ok>')% File exists.
else
     fprintf('KPOINTS<error>')% File does not exist.
end
if isfile('EIGENVAL')
     fprintf('EIGENVAL<ok>')% File exists.
else
     fprintf('EIGENVAL<error>')% File does not exist.
end
if isfile('DOSCAR')
     fprintf('DOSCAR<ok>')% File exists.
else
     fprintf('DOSCAR<error>')% File does not exist.
end
%% for a default print,just type 
tmpfig = bandplot;
%% for a user-defined bandplot
EIGENCAR = EIGENVAL_read();
bandplot(EIGENCAR,[-6,6],'title','Band Stucture CuI w SOC','Color','b');

