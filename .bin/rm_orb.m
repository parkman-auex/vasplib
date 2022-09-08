%% Remove orbitals of H_xyz

function [Out_H_xyz,new_orb] = rm_orb(H_xyz,orb_rm_index_list,orb_list,mode)


% x = linspace(1,120,120);          % orbital number
% y = [51,52,55,56,111,112,115,116]; % select pxpy orbital
% z = setdiff(x,y);   % the orbitals to delet
%
if exist('POSCAR','file')
    %     Ns = [1 0 0;0 1 0;0 0 1];
    %     Ns(fin_dir,:) = Ns(fin_dir,:) * repeatnum;
    POSCAR_read;
    %     fin_dir_list = [0 0 0];
    %     if strcmp(glue_edges,'false')
    %         fin_dir_list(fin_dir) = 1;
    %     end
    %     [Rm,~] = supercell(Ns,Rm,sites,Atom_name,Atom_num,fin_dir_list,'POSCAR_super_fin');
end

%
if nargin < 4
    mode = 't';
end
%% init
if strcmp(mode,'t')
    [NRPTS,~]=size(H_xyz);
    WAN_NUM=length(H_xyz(1).Hnum);
    Out_H_xyz = H_xyz;
elseif strcmp(mode,'c')
    Hnum_list = H_xyz.HnumL ;
    vector_list = H_xyz.vectorL ;
    [NRPTS,~]=size(vector_list);
    %    WAN_NUM=length(Hnum_list(:,:,1));
    WAN_NUM_2 = length(Hnum_list(:,1));
    WAN_NUM=sqrt(WAN_NUM_2);
elseif strcmp(mode,'s')
    Hnum_list = H_xyz.HnumL ;
    vector_list = H_xyz.vectorL ;
    [NRPTS,~]=size(vector_list);
    WAN_NUM=length(Hnum_list{1});
end

%% nargin
if nargin < 3
    orb_list =zeros(WAN_NUM,3);
end
rm_orb_num = length(orb_rm_index_list);
fprintf('Total num of rm_orb:%d', rm_orb_num);
new_orb  = orb_list;
new_orb(orb_rm_index_list,:) = [];

if strcmp(mode,'t')
    %disp(orb_rm_index_list);
    %         disp(length(new_orb));
    %                    disp(new_orb);
    %% By cchen
    for k = 1:NRPTS   % number of home cells
        %     x=[2,6,10,14,18,22,26,30,34,38,42,46];  %
        Out_H_xyz(k).Hnum(orb_rm_index_list,:)=[];
        Out_H_xyz(k).Hnum(:,orb_rm_index_list)=[];
        Out_H_xyz(k).Hcoe(orb_rm_index_list,:)=[];
        Out_H_xyz(k).Hcoe(:,orb_rm_index_list)=[];
        
        %     Out_H_xyz(k).Hnum(orb_rm_index_list,:)=[];
        %     Out_H_xyz(k).Hnum(:,orb_rm_index_list)=[];
        %     Out_H_xyz(k).Hcoe(orb_rm_index_list,:)=[];
        %     Out_H_xyz(k).Hcoe(:,orb_rm_index_list)=[];
    end
elseif strcmp(mode,'s')
    %disp(orb_rm_index_list);
    %         disp(length(new_orb));
    %                    disp(new_orb);
        disp('sparse_mode');
    for k = 1:NRPTS   % number of home cells
        OUT_WAN_NUM  = WAN_NUM-rm_orb_num ;
        Hnum_list{k}(orb_rm_index_list,:)=[];
        Hnum_list{k}(:,orb_rm_index_list)=[];
        OUT_Hnum_list = Hnum_list;
        OUT_vector_list = vector_list ;
        Out_H_xyz.vectorL  = OUT_vector_list ;
        Out_H_xyz.HnumL = OUT_Hnum_list;
    end
    
elseif strcmp(mode,'c')
    disp('concise_mode');
    OUT_WAN_NUM  = WAN_NUM-rm_orb_num ;
    OUT_WAN_NUM_2 = OUT_WAN_NUM *2;
    orb_rm_index_list_concise = concise_rm_orb_list(orb_rm_index_list,WAN_NUM);
    Hnum_list( orb_rm_index_list_concise ,:)=[];
    OUT_Hnum_list = Hnum_list;
    OUT_vector_list = vector_list ;
    Out_H_xyz.vectorL  = OUT_vector_list ;
    Out_H_xyz.HnumL = OUT_Hnum_list;
elseif strcmp(mode,'sym')
    disp('sym_mode');
end


end

function orb_rm_index_list_concise = concise_rm_orb_list(orb_rm_index_list,WAN_NUM)
orb_rm_index_list_concise = [];
rm_orb_num = length(orb_rm_index_list);
for i=1:rm_orb_num
    temp_list1 = (orb_rm_index_list(i)-1)*WAN_NUM+[1:WAN_NUM];
    temp_list2 = orb_rm_index_list(i)+[0:WAN_NUM-1]*WAN_NUM;
    orb_rm_index_list_concise  = [orb_rm_index_list_concise ,temp_list1,temp_list2 ];
end
end