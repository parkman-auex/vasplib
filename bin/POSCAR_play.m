%% note: POSCAR is Direct mode
% input : POSCAR in current dir
% Output : cystal_information: 
%                             a_crystal_constance
%                             a1,a2,a3,Rm=[a1;a2;a3]
%                             atom_name;atom_num
%                             coordinates pattern
%          element_information:
%                             sites 
% note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
% Rm
    %a1 a1x a1y a1z
    %a2 a2x a2y a2z
    %a3 a3x a3y a3z
% unit A?e-6
function POSCAR_play(sites)
sites_num = length(sites);
    for i=1:sites_num
            fprintf("%d ",i);
            fprintf("%4s ",sites(i).name);
            fprintf("%f ",sites(i).rc1);
            fprintf("%f ",sites(i).rc2);
            fprintf("%f ",sites(i).rc3);
            fprintf("\n")

    end
end
    
    
