%% POSCAR_readin
%
% get poscar information from vasp and others
% * Label:
%
%% Description of the Function:
%%
%% Usage:
%
% * [Rm,sites,Atom_name,Atom_num]= POSCAR_readin(filename,mode)
% * [Rm,sites,Atom_name,Atom_num] = POSCAR_readin(filename)
% * [Rm,sites,Atom_name,Atom_num] = POSCAR_readin()
%
%% Input:
%
% # POSCAR in current dir
%
%% Output:
%
% # a_crystal_constance
% # Rm
% # Atom_name
% # Atom_num
% # element_information: sites
% unit A?e-6
% ## seq: a local label of an ion
% ## inseq: nth ion for the elements
% ## rc1 rc2 rc3 : Direct fractional coordinates
% ## nameseq :the local label of the element
% ## name : a string to
% ## mag:  the local mag
% ## ion_type: waiting, need a tanslator
% ## orbs: l quantum_number or 's' 'p' 'd'
% ## orb_sym: a batch of function {'1','x','y','z'} d & f is waiting
%% example:
%   commmad
%
%
%   result
%
%% Note:
%
%  site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[],'mag',[],'ion_type',[],'orb',[],'orb_sym',[]);
%
%% Change log
%
% * Document Date: 2020/12/03
% * Creation Date: 2020/12/03
% * Last updated : 2020/12/03
% * Last updated : 2020/12/25
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code :
%
%% note: POSCAR is Direct mode
% input : POSCAR in current dir


function [Rm,sites,Atom_name,Atom_num,elements,a_crystal_constance]=POSCAR_readin(filename,mode,options)
%--------  narg  --------
arguments
    filename = 'POSCAR';
    mode = 'vasp';
    options.digits = 6;
end
warning("This function is deprecated and will be removed in the future, please use POSCAR_read instead")
return



%--------  init  --------
elements = readtable('elements.txt');
%--------  init  --------
digits(options.digits); %take care
%--------  chek  --------
if ~exist(filename,'file')
    fprintf('No such file: %s\n',filename);
    error('file not exist!');
end
%--------  fbug  --------
% disp(mode);

%--------  juge  --------
switch mode
    case 'vasp'
        site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
        formatSpec = '%s%s%s%s%s%[^\n\r]';
    case 'tbsk'
        site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[],'orb',[],'orb_sym',[]);
        formatSpec = '%s%s%s%s%s%s%[^\n\r]';
    case 'tbsym'
        site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[],'element',[],'orb',[],'orb_sym',[],'spin',[]);
        formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';
    case 'list'
        site=struct('seq',[],'rc1',[],'rc2',[],'rc3',[],'Hue',[],'surf_level',[],'hing_level',[]);
    
end

POSCAR = POSCAR_cell_read(filename,formatSpec);
a_crystal_constance=str2double(char(POSCAR(1,1)));
% Rm
a1=[str2double(char(POSCAR(2,1))) str2double(char(POSCAR(2,2))) str2double(char(POSCAR(2,3)))];
a2=[str2double(char(POSCAR(3,1))) str2double(char(POSCAR(3,2))) str2double(char(POSCAR(3,3)))];
a3=[str2double(char(POSCAR(4,1))) str2double(char(POSCAR(4,2))) str2double(char(POSCAR(4,3)))];
Rm=[a1;a2;a3];
% read atom
n_atom=length(POSCAR(5,:));
for i=1:n_atom
    if POSCAR(5,i) ~= ""
        Atom_name(i) = POSCAR(5,i);
    end
end
% need
for i=1:length(Atom_name)
    Atom_num(i)=str2double(char(POSCAR(6,i)));
end
%site_num
sites_num=sum(Atom_num);
sites=repmat(site,[1 sites_num]);
% coordinate pattern
% Coordinates_pattern=POSCAR(7,1);
% temp : frac support
% struct
sequence=1;
n_atom=length(Atom_name);
labelcut_list = labelcut_list_gen(Atom_num);

%first
if n_atom >= 1
    for i=1:n_atom
        inseq =1;
        for j=labelcut_list(i,1):labelcut_list(i,2)
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
            sites(sequence).rc1=str2double(char(POSCAR(j,1)));
            sites(sequence).rc2=str2double(char(POSCAR(j,2)));
            sites(sequence).rc3=str2double(char(POSCAR(j,3)));
            if strcmp(mode,'tbsk')
                sites(sequence).orb = string(POSCAR(j,5));
                sites(sequence).orb_sym = sym(POSCAR(j,6));
            elseif strcmp(mode,'tbsym')
                sites(sequence).element = string(POSCAR(j,4));
                sites(sequence).orb = string(POSCAR(j,5));
                sites(sequence).orb_sym = string(POSCAR(j,6));
                sites(sequence).spin = str2double(POSCAR(j,7));
            end
            %
            sequence=sequence+1;
            inseq =inseq + 1;
        end
    end
end

if length(POSCAR)>j & strcmp(mode,'mag')
    %disp('mag_mode');
    sequence=1;
    beginline=j+1;
    %first
    for j=beginline:beginline+Atom_num(1)-1
        %id
        sites(sequence).mag=double(POSCAR(j,1:3));
        sequence=sequence+1;
    end
    %other
    if n_atom >= 2
        for i=2:n_atom
            beginline=beginline+Atom_num(i-1);
            for j=beginline:beginline+Atom_num(i)-1
                %id
                sites(sequence).mag=double(POSCAR(j,1:3));
                sequence=sequence+1;
            end
        end
    end
end


end

function labelcut_list = labelcut_list_gen(Atom_num)
    n = length(Atom_num);
    beginline =8 ;
    sum_n =beginline;
    for i =1:n
        sum_n2 = sum_n+Atom_num(i)-1;
        labelcut_list(i,:) = [sum_n sum_n2];
        sum_n = sum_n2+1;
    end

end
function POSCAR = POSCAR_cell_read(filename,formatSpec)
    startRow = 2;
    delimiter = {'\t',' '};
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    rawStringColumns = string(raw(:,:));
    POSCAR=rawStringColumns;
    clearvars filename delimiter startRow formatSpec fileID dataArray ans rawStringColumns raw;
end

