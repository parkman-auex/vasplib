function [Rm_standard,sites_standard,Atom_name,Atom_num_standard] = standardize_cell(POSCAR_name,symprec,to_primitive,no_idealize,spglib_path,spglib_include)
if nargin < 1
    POSCAR_name = 'POSCAR';
end
if nargin < 2 
    symprec = 0.01;
end
if nargin < 3
    to_primitive = 0 ;
end
if nargin < 4
    no_idealize  = 0 ;
end
if nargin < 5
    spglib_path = '/usr/local/lib/';
end
if nargin < 6
    spglib_include = '/usr/local/include/';
end


warning off;
if not(libisloaded('libsymspg'))
    addpath(spglib_path);
    addpath(spglib_include);
    loadlibrary('libsymspg','spglib.h');
end

[Rm,sites,Atom_name,Atom_num] = POSCAR_read(POSCAR_name);

lattice = Rm;
num_atom = sum(Atom_num);
position = [[sites.rc1].',[sites.rc2].',[sites.rc3].'].';

types = [];
for i = 1:length(Atom_num)
    types = [types;ones(Atom_num(i),1)*i];
end

[num_atom,Rm_standard,position_standard,types_standard]=calllib('libsymspg','spg_standardize_cell',...
                     lattice,...
                     position,...
                     types,...
                     num_atom,...
                     to_primitive,...
                     no_idealize,...
                     symprec);
position_standard(:,num_atom+1:end) = [];
types_standard(num_atom+1:end) = [];
sites_standard=sites(1:num_atom);
Atom_num_standard = zeros(size(Atom_name,2),1);
for i =1:num_atom
    sites_standard(i).nameseq = types_standard(i);
    sites_standard(i).rc1 = position_standard(1,i);
    sites_standard(i).rc2 = position_standard(2,i);
    sites_standard(i).rc3 = position_standard(3,i);
    Atom_num_standard(types_standard(i)) = Atom_num_standard(types_standard(i)) +1;
    sites_standard(i).inseq = Atom_num_standard(types_standard(i));
    sites_standard(i).name = strcat(Atom_name(types_standard(i)),string(sites_standard(i).inseq));
end
POSCAR_gen(Rm_standard,sites_standard,Atom_name,Atom_num_standard,'PPOSCAR');

end

