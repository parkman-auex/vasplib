%% make sure present work and save the var
% usage: dirstring=init_v2w(level_cut,search_range,mesh)
%        dirstring=init_v2w(level_cut)
%        dirstring=init_v2w()
%% Attention 
%  input :
%         INCAR;POSCAR;KPOINT(band);EIGENVAL(band);DOSCAR(band);
%          SOCflag;
%clear;

function dirstring=init_v2w(SOCFLAG,Erange,options)
arguments
    SOCFLAG = 0;
    Erange = [-3,3];
    options.check = false;
    options.Projector_list = [];
end
import park.*;
SOCflag=SOCFLAG;
[Rm,sites,Atom_name,Atom_num]=POSCAR_readin();
try
    check = true;
    EIGENCAR=EIGENVAL_read();
catch
    EIGENCAR = [];
    check = false;
end
if options.check && nargin < 2
%% while
    while(~strcmp(input('y/n for continue or break to findbands','s'),'n'))
       findbands(EIGENCAR); 
    end
elseif options.check
    findbands(EIGENCAR,Erange);
else
end
%
kpath_card_gen();
%
if isempty(options.Projector_list)
    num_wan=write_pj();
else
    num_wan=write_pj(options.Projector_list);
end

wt_in_gen(Rm,sites,Atom_num,Atom_name,SOCflag);
if check && options.check
    [num_bands,Nmin,Nmax] = selectbands(EIGENCAR,Rm);
    Nbands = size(EIGENCAR,1);
    wannier90_win_gen(SOCflag,num_bands,num_wan,Nmin,Nmax,Nbands);
else
    wannier90_win_gen(SOCflag,10,num_wan,1,10,10);
end

disp("all done");

end
