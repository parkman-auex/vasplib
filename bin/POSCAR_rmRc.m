%% Remove a specfic rc or a rc list in POSCAR
% usage : [Rm,sites] = POSCAR_rmRc(rc_rm_list,Accuracy,filename,Rm,sites,Atom_name,Atom_num)
function [Rm,sites] = POSCAR_rmRc(rc_rm_list,Accuracy,filename,Rm,sites,Atom_name,Atom_num)


%function [Rm,sites] = POSCAR_gen(Rm,sites,Atom_name,Atom_num,filename)

%% nargin
if nargin < 4
    POSCAR_read;
end
if nargin < 3
    filename = 'POSCAR_rmRc';
end
if nargin < 2
   Accuracy = 4;
end


Rc_list = [[sites.rc1]',[sites.rc2]',[sites.rc3]'];

%% same Accuracy
rc_rm_list = round(rc_rm_list.*10^Accuracy)./10^Accuracy ; 
% disp(rc_rm_list);
rc_rm_list = unique(rc_rm_list,'rows');
% disp(rc_rm_list);

Rc_list  = round(Rc_list.*10^Accuracy)./10^Accuracy;

seq_list = find_unique_list_in_a_list(rc_rm_list,Rc_list);
% disp(seq_list);
% check each type atom losenum 
atomtype_list = [sites(seq_list).nameseq]';
%disp(atomtype_list );
% creat atom_type list
Atom_type_list = [1:length(Atom_name)]';
%disp(Atom_type_list)
atomtype_list = [Atom_type_list;atomtype_list];
%[all_one,seq_list_init]=ismember(atomtype_list,Atom_type_list','rows');



repeat_list= histc(atomtype_list,unique(atomtype_list))-1;
%disp(repeat_list );
% delete atom number
Atom_num = Atom_num - repeat_list';


% findif0
%disp(Atom_num);
[atom_unique,atom_seq_list]=ismember([0],Atom_num' ,'rows');
if sum(atom_unique) > 0
Atom_num(atom_seq_list) = [];
Atom_name(atom_seq_list) = [];
end

%% rm rc
sites(seq_list) = [];

[Rm,sites] = POSCAR_gen(Rm,sites,Atom_name,Atom_num,filename);
end

function seq_list = find_unique_list_in_a_list(list_a,list_b)
    seq_list = [];
    %listb>lista
    if length(list_b) <length(list_a)
        warning('list_2nd longer than list_1st')
    end
    if isempty(list_a)
        disp('No rm list')
        return
    end
    [all_one,seq_list_init]=ismember(list_a,list_b,'rows');
    %unique
%     [seq_unique,seq_list_init_unique] = unique(seq_list_init,'rows');
%     all_one = all_one(seq_list_init_unique);
    nseq_list = length(all_one);

    for i =1:nseq_list
        if all_one(i) == 1
            seq_list= [seq_list;seq_list_init(i,:)];
        end
    end


end