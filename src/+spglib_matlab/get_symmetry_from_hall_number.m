function [Nsym,rotation,translation] = get_symmetry_from_hall_number(hall_number,spglib_path,spglib_include)
if nargin < 2
    spglib_path = '/usr/local/lib/';
end
if nargin < 3
    spglib_include = '/usr/local/include/';
end
 max_size = 192;
 rotation =  libpointer('int32Ptr',zeros(3,3,max_size ));
 translation = libpointer('doublePtr',zeros(3,max_size ));
%[hall_number, rotation ,  translation]  =
 [Nsym, rotation ,  translation]  = calllib('libsymspg','spg_get_symmetry_from_database',rotation, translation,hall_number);
                  rotation = reshape(rotation,3,3,[]);
                 for i = 1:Nsym
                     rotation(:,:,i) =rotation(:,:,i).';
                 end
                 rotation(:,:,Nsym+1:end) = [] ;
                 translation = translation.';
                 translation(Nsym+1:end,:) = [] ;
end

