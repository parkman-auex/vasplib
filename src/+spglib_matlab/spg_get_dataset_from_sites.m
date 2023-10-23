function SpglibDataset  = spg_get_dataset_from_sites(lattice,position,types,symprec,spglib_path,spglib_include)
if nargin < 3
    POSCAR_name = 'POSCAR';
    [Rm,sites,~,Atom_num] = POSCAR_readin(POSCAR_name);
    lattice = Rm;
    position = [[sites.rc1].',[sites.rc2].',[sites.rc3].'].';
    types = [];
    for i = 1:length(Atom_num)
        types = [types,ones(1,Atom_num(i)*i)];
    end
end
if nargin < 4 
    symprec = 0.01;
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
num_atom = length(types);



% SpglibDatasetPtr = libpointer();
[SpglibDatasetPtr,~,~,~] = calllib('libsymspg','spg_get_dataset',lattice,position,types,num_atom,symprec);
%SpglibDataset = calllib('libsymspg','spg_free_dataset',SpglibDatasetPtr);
SpglibDataset = SpglibDatasetPtr.Value;
clear SpglibDatasetPtr;
%% ptr
% rotations
rotationsPtr = SpglibDataset.error0;
rotationsPtr.setdatatype('int32Ptr',3,3*SpglibDataset.n_operations);
rotations = reshape(rotationsPtr.Value,3,3,[]);
                 for i = 1:SpglibDataset.n_operations
                     rotations(:,:,i) =rotations(:,:,i).';
                 end
SpglibDataset.rotations  =    rotations ;
% transformation
translationsPtr = SpglibDataset.error1;
translationsPtr.setdatatype('doublePtr',3,SpglibDataset.n_operations);
translations = translationsPtr.Value';
SpglibDataset.translations  =    reshape(translations,1,3,SpglibDataset.n_operations);
% wyckoffs
wyckoffsPtr = SpglibDataset.wyckoffs;
wyckoffsPtr.setdatatype(wyckoffsPtr.DataType,SpglibDataset.n_atoms)
SpglibDataset.wyckoffs = wyckoffsPtr.Value;
% site_symmetry_symbols
site_symmetry_symbolsPtr = SpglibDataset.error2;
site_symmetry_symbolsPtr.setdatatype('int8Ptr',7,SpglibDataset.n_atoms);
SpglibDataset.site_symmetry_symbols = string(char(site_symmetry_symbolsPtr.Value.'));
% equivalent_atoms
equivalent_atomsPtr = SpglibDataset.equivalent_atoms;
equivalent_atomsPtr.setdatatype(equivalent_atomsPtr.DataType,SpglibDataset.n_atoms)
SpglibDataset.equivalent_atoms = equivalent_atomsPtr.Value;
% crystallographic_orbits
crystallographic_orbitsPtr = SpglibDataset.crystallographic_orbits;
crystallographic_orbitsPtr.setdatatype(crystallographic_orbitsPtr.DataType,SpglibDataset.n_atoms)
SpglibDataset.crystallographic_orbits = crystallographic_orbitsPtr.Value;
% std_types
std_typesPtr = SpglibDataset.std_types;
std_typesPtr.setdatatype(std_typesPtr.DataType,SpglibDataset.n_atoms)
SpglibDataset.std_types = std_typesPtr.Value;
% std_positions
std_positionsPtr = SpglibDataset.error3;
std_positionsPtr.setdatatype('doublePtr',3,SpglibDataset.n_atoms);
std_positions =std_positionsPtr.Value';
SpglibDataset.std_positions  =    std_positions;
% point_group_symbole
SpglibDataset.pointgroup_symbol = string(strcat(char(SpglibDataset.pointgroup_symbol)));
SpglibDataset.transformation_matrix = int32(reshape(SpglibDataset.transformation_matrix,3,3));
SpglibDataset.primitive_lattice =(reshape(SpglibDataset.primitive_lattice,3,3));
SpglibDataset.std_lattice =(reshape(SpglibDataset.std_lattice,3,3));
SpglibDataset.std_rotation_matrix =(reshape(SpglibDataset.std_rotation_matrix,3,3));
% 
SpglibDataset = rmfield(SpglibDataset,{'error0','error1','error2','error3','mapping_to_primitive','std_mapping_to_primitive'});
SpglibDataset.international_symbol = string(strcat(char(SpglibDataset.international_symbol)));
SpglibDataset.hall_symbol = string(strcat(char(SpglibDataset.hall_symbol)));
end

% function
% end