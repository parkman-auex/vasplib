%% make sure present work and save the var
% usage: dirstring=vasp_matlab_prework(level_cut,search_range,mesh)
%        dirstring=vasp_matlab_prework(level_cut)
%        dirstring=vasp_matlab_prework()

%clear;
%% make sure , we have 
function dirstring=vasp_matlab_prework(level_cut,search_range,mesh)
if nargin==2
    mesh=[13 13 13];
end
if nargin==1
    mesh=[13 13 13];
    search_range=[1 1 1];
end
if nargin==0
    level_cut=3;
    mesh=[13 13 13];
    search_range=[1 1 1];
end

POSCAR_read;
[Atom_store,nn_store,Rnn]=nn(Rm,sites,search_range);
H=H_TB_gen(level_cut,nn_store,sites);
[klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm);
%

[klist_dos]=kmesh3D(Rm,mesh);

% save 
    dirname=pwd;
    dirname=strsplit(dirname,'/');
    dirstring=dirname(length(dirname))+".mat";
    save(char(dirstring));
end

