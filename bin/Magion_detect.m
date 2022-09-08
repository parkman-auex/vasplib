function [mag_ion_seq,sites,Atom_name,Atom_num] = Magion_detect(mag_element_seq_list,Rm,sites,Atom_name,Atom_num,elements)
%UNTITLED2 此处显示有关此函数的摘要

% nargin
if nargin <1
    mag_element_seq_list = 56:72; % rare_earth add more!!
end
if nargin <2
    % read POSCAR
    [Rm,sites,Atom_name,Atom_num,elements]=POSCAR_readin('POSCAR','vasp');
end
% make peridoic table
elements.Properties.RowNames = elements.atom_symbol;
%
fprintf('We have %d ions:\n',sum(Atom_num));
mag_ion_seq  = -1;

for i = 1:length(Atom_name)
    tmpdata = table2array(elements(Atom_name(i),{'atom_number','n'}));
    element_seq = tmpdata(1);
    fprintf('    %d %s(%d) ions \n',Atom_num(i),Atom_name(i),element_seq);
    if ismember(element_seq,mag_element_seq_list)
        if mag_ion_seq >0
            error('Dont support multiple rare-earth element at present');
        end
        if Atom_num(i) > 1
            warning('Dont support more rare-earth ions in unit cell  at present, there may have some bugs');
        end
        mag_ion_seq = i;
        mag_ion_site_seq = sum(Atom_num(1:mag_ion_seq));
    end
end
if mag_ion_seq == -1
    fprintf('do not detect re-earth element');
    error('!!');
else
    fprintf('Rare earth element: %s\n',Atom_name(mag_ion_seq));
end
% distri
Atom_name_mag = Atom_name(mag_ion_seq);
Atom_num_mag  = Atom_num(mag_ion_seq);
sites_mag     = sites(mag_ion_site_seq);
POSCAR_gen(Rm,sites_mag,Atom_name_mag,Atom_num_mag,'POSCAR_mag');
Atom_name(mag_ion_seq) = [];
Atom_num(mag_ion_seq) = [];
sites(mag_ion_site_seq) = [];
POSCAR_gen(Rm,sites,Atom_name,Atom_num,'POSCAR_others');

sites = [sites_mag,sites];
Atom_name = [Atom_name_mag,Atom_name];
Atom_num = [Atom_num_mag,Atom_num];
%   此处显示详细说明
end

