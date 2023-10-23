function [SG_number,SG_international] = get_international(POSCAR_name,symprec,spglib_path,spglib_include)
if nargin < 1
    POSCAR_name = 'POSCAR';
end
if nargin < 2 
    symprec = 0.01;
end
if nargin < 3
    spglib_path = '/usr/local/lib/';
end
if nargin < 4
    spglib_include = '/usr/local/include/';
end


warning off;
if not(libisloaded('libsymspg'))
    addpath(spglib_path);
    addpath(spglib_include);
    loadlibrary('libsymspg','spglib.h');
end

[Rm,sites,~,Atom_num] = POSCAR_readin(POSCAR_name);

lattice = Rm;
position = [[sites.rc1].',[sites.rc2].',[sites.rc3].'].';
types = [];
for i = 1:length(Atom_num)
    types = [types;ones(1,Atom_num(i)*i)];
end
num_atom = sum(Atom_num);
SG_international = libpointer('int8Ptr',zeros(11,1));
[SG_number,SG_international,~,~,~]=calllib('libsymspg','spg_get_international',SG_international,...
                     lattice,...
                     position,...
                     types,...
                     num_atom,...
                     symprec);
SG_international = strcat(char(SG_international).');

fprintf('space_group_type: %s\n',SG_international);
fprintf('space_group_number: %d\n',SG_number);
end

