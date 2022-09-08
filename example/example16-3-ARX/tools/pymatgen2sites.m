function [Rm,sites,Atom_name,Atom_num] = pymatgen2sites(pymatgen_struct,filename)
if nargin <2
    mode = 'nogen';
else
    mode = 'gen';
end
    lattice = pymatgen_struct.lattice;
    sites_temp = pymatgen_struct.sites;
    Rm = lattice.matrix;
    NIONS = length(sites_temp);
    for i = 1:NIONS
        Atom_list{i} = sites_temp(i).label;
    end
    [Atom_list_cell,cutlabel_list] = unique(Atom_list);
    [cutlabel_list,reseq] = sort(cutlabel_list);
    NELEMENTS = length(Atom_list_cell );
    %cutlabel_list = cutlabel_list(reseq);
    cutlabel_list(NELEMENTS+1) = NIONS+1;
    for i =1:NELEMENTS 
        Atom_name(i) = string(Atom_list_cell{i});
        Atom_num(i) =  cutlabel_list(i+1)-cutlabel_list(i);
    end
    Atom_name = Atom_name(reseq);
    site = struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
    sites=repmat(site,[1 NIONS]);
    sequence=1;
    if NELEMENTS >= 1
        for i=1:NELEMENTS
            inseq =1;
            for j= cutlabel_list(i):cutlabel_list(i+1)-1
                %id
                sites(sequence).seq=sequence;
                %inid
                sites(sequence).inseq=inseq ;
                %name
                sites(sequence).nameseq=i;
                sites(sequence).name=Atom_name(i)+num2str(sites(sequence).inseq);
                sites(sequence).ion_num=Atom_num(i);
                %sites(sequence).ion_type_num=i;
                %incoordinate
                sites(sequence).rc1=sites_temp(j).abc(1);
                sites(sequence).rc2=sites_temp(j).abc(2);
                sites(sequence).rc3=sites_temp(j).abc(3);
                %
                sequence=sequence+1;
                inseq =inseq + 1;
            end
        end
    end
    if strcmp(mode,'gen')
         POSCAR_gen(Rm,sites,Atom_name,Atom_num,filename);
    end
end