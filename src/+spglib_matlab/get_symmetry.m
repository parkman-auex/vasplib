function [Nsym,rotation,translation] = get_symmetry(POSCAR_name,symprec,spglib_path,spglib_include)
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

[Rm,sites,~,Atom_num] = POSCAR_read(POSCAR_name);
max_size = 192;
rotation =  libpointer('int32Ptr',zeros(3,3,max_size));
translation = libpointer('doublePtr',zeros(3,max_size));
lattice = Rm;
position = [[sites.rc1].',[sites.rc2].',[sites.rc3].'].';
types = [];
for i = 1:length(Atom_num)
    types = [types;ones(1,Atom_num(i)*i)];
end
num_atom = sum(Atom_num);
[Nsym,rotation,translation,lattice,position,types] = calllib('libsymspg','spg_get_symmetry',rotation,...
                     translation,...
                     max_size,...
                     lattice,...
                     position,...
                     types,...
                     num_atom,...
                     symprec);
                 
                 rotation = reshape(rotation,3,3,[]);
                 for i = 1:Nsym
                     rotation(:,:,i) =rotation(:,:,i).';
                 end
                 rotation(:,:,Nsym+1:end) = [] ;
                 translation = translation.';
                 translation(Nsym+1:end,:) = [] ;
fprintf('# %3d space_group_operations:\n',Nsym);
for i = 1:Nsym
    fprintf('* rotation: # %d\n',i)
    fprintf('%3d %3d %3d\n%3d %3d %3d\n%3d %3d %3d\n',...
        rotation(1,1,i),rotation(1,2,i),rotation(1,3,i),...
        rotation(2,1,i),rotation(2,2,i),rotation(2,3,i),...
        rotation(3,1,i),rotation(3,2,i),rotation(3,3,i)...
        );
    fprintf('  ---------\n');
    fprintf('  %3.1f %3.1f %3.1f\n',translation(i,1),translation(i,2),translation(i,3));
end
end

