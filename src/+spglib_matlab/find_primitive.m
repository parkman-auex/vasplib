function [Rm_primitive,sites_primitive,Atom_name,Atom_num_primitive] = find_primitive(POSCAR_name,symprec,spglib_path,spglib_include)
to_primitive = 1;
no_idealize = 0;
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
[Rm_primitive,sites_primitive,Atom_name,Atom_num_primitive] = standardize_cell(POSCAR_name,symprec,to_primitive,no_idealize,spglib_path,spglib_include);
end