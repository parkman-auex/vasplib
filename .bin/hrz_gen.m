%% Hr_sprse2H_hrz
%
% change a Hr_sparse into H_hrz by a certain Kpoints and direction
% useful for cut_piece and GreenFunction
% 
%
% * Label: hr
%
%% Description of the Function:
%%
%% Usage: 
%
% * [H_hrz,Hnum_list_hrz,vector_list_hrz]= hrz_gen(hr_sparse,kpoints_f,fin_dir)
% * H_hrz = hrz_gen(hr_sparse,kpoints_f,fin_dir)
% 
%% Input:
%  
% # input1:
% # input2:
% # input3:
%
%% Output:
%
% # fig:
% # ax:
%
%% example:
%
%   
%% Note: 
%
%  Take advantage of the scope of application of the function.
%
%% Change log
%
% * Document Date: 2020/12/16
% * Creation Date: 2020/12/16
% * Last updated : 2020/12/16
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code : 
%
function [H_hrz,Hnum_list_hrz,vector_list_hrz]= hrz_gen(hr_sparse,kpoints_f,fin_dir,mode)
%--------  narg  --------
    if nargin < 4
       mode = '1D';
    end

%--------  init  --------
    import vasplib_tool.*
    Hnum_list = hr_sparse.HnumL;
    WAN_NUM=length(Hnum_list{1});
    NRPTS=length(Hnum_list);
    vector_list = hr_sparse.vectorL;
    
    factor_list = exp(1i*2*pi*(vector_list*kpoints_f.'));
    Hnum_list3 = Hnum_list3_gen(Hnum_list);
    Hnum_kf_list(:,:,NRPTS) = Hnum_list3(:,:,NRPTS);
    for i = 1 : NRPTS
        Hnum_kf_list(:,:,i) = Hnum_list3(:,:,i)*factor_list(i);
    end
%--------  main  --------
switch mode
    case '0D'
        H_hrz = sum(Hnum_kf_list,3);
        Hnum_list_hrz = Hnum_kf_list ;
        vector_list_hrz = vector_list;
    case '1D'
        [vector_list_new,sort_label] = sortrows(vector_list,fin_dir) ;% sort fin_dir
        Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label); % sort hr
        [unique_dir,unique_label]= unique(vector_list_new(:,fin_dir));
        %cutlist
        cutlist(:,1)= unique_label;
        cutlist(1:end-1,2)= unique_label(2:end)-1;
        cutlist(end,2) = NRPTS;
        %--------  H_hrz  --------
        vector_list_hrz = zeros(length(unique_dir),3);
        NRPTS_hrz = length(unique_label);
        Hnum_list_hrz = zeros(WAN_NUM,WAN_NUM,NRPTS_hrz);
        for i = 1:NRPTS_hrz
            vector_list_hrz(i,fin_dir) = unique_dir(i,:);
            Hnum_list_hrz(:,:,i) = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
        end
        H_hrz.HnumL = Hnum_list_hrz;
        H_hrz.vectorL = vector_list_hrz;
    case '2D'
         [vector_list_new,sort_label] = sortrows(vector_list,fin_dir) ;% sort fin_dir
         [unique_dir,unique_label]= unique(vector_list_new(:,fin_dir),'rows');
         Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label); % sort hr
         %cutlist
         cutlist(:,1)= unique_label;
         cutlist(1:end-1,2)= unique_label(2:end)-1;
         cutlist(end,2) = NRPTS;
         %--------  H_hrz  --------
         vector_list_hrz = zeros(length(unique_dir),3);
         NRPTS_hrz = length(unique_label);
         Hnum_list_hrz = zeros(WAN_NUM,WAN_NUM,NRPTS_hrz);
         for i = 1:NRPTS_hrz
            vector_list_hrz(i,fin_dir) = unique_dir(i,:);
            Hnum_list_hrz(:,:,i) = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
         end
         H_hrz.HnumL = Hnum_list_hrz;
         H_hrz.vectorL = vector_list_hrz;
end
end

function H = Hnum_list3_gen(Hnum_list)
    WAN_NUM=length(Hnum_list{1});
    NRPTS=length(Hnum_list);
    H = reshape(full(cell2mat(Hnum_list)),WAN_NUM,WAN_NUM,NRPTS);
end