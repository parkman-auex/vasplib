function H_E = add_electric_field(E)
arguments
    E (1,3) {mustBeReal} = [1 0 0] % V/Ang, [Ex Ey Ez]
end
try
    xyz = readtable('wannier90_centres.xyz','NumHeaderLines',2,'FileType','text');
catch
    error('Can not read wannier90_centres.xyz');
end
wcentres = xyz( strcmp(xyz.Var1,{'X'}), 2:4);
E_onsite = table2array(wcentres) * E';
H_E = diag(E_onsite);
end