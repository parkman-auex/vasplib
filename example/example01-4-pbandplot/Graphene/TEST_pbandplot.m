%% make sure you have copy the POSCAR KPOINTS EIGENVAL DOSCAR in this directory
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
if isfile('PBAND_C.dat')
     fprintf('PBAND_C.dat<ok>')% File exists.
else
     fprintf('PBAND_C.dat<error>')% File does not exist.
end
%% for a default print,just type 
vasplib.pbandplot('Ecut',[-6,6]);

