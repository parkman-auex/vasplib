function H_hr = hr_converter(mode,write)
% see function 'add_zeeman_field' for more infomation
arguments
    mode {mustBeMember(mode,{'old_to_new','new_to_old'})} = 'old_to_new';
    write logical = true;
end
H_hr = HR.from_wannier90();
HnumL_tmp = zeros(size(H_hr.HnumL));
Norb = H_hr.Nbands/2;
if strcmp(mode , 'old_to_new')
    for i = 1:Norb
        for j = 1:Norb
            HnumL_tmp(2*i-1,2*j-1,:) = H_hr.HnumL(i     ,j     ,:);
            HnumL_tmp(2*i  ,2*j-1,:) = H_hr.HnumL(i+Norb,j     ,:);
            HnumL_tmp(2*i-1,2*j  ,:) = H_hr.HnumL(i     ,j+Norb,:);
            HnumL_tmp(2*i  ,2*j  ,:) = H_hr.HnumL(i+Norb,j+Norb,:);
        end
    end
elseif strcmp(mode , 'new_to_old')
    for i = 1:Norb
        for j = 1:Norb
            HnumL_tmp(i     ,j     ,:) = H_hr.HnumL(2*i-1,2*j-1,:);
            HnumL_tmp(i+Norb,j+Norb,:) = H_hr.HnumL(2*i  ,2*j  ,:);
        end
    end
end
H_hr.HnumL = HnumL_tmp;
if write
    H_hr.Gen_hr('wannier90_hr.converted.dat');
end
end