%% hr2hr_sparse
%% use for construct sparse matrics hr, 
% suppose all Degen = 1
% with a static struct 
% usage : [hr_concise,Hnum_list,vector_list] = hr2hr_sparse(H_xyz,mode)
% we will define a 3D mat a lisr_for NRPT list(V) and a list for Degen()

function [hr_sparse,Hnum_list,vector_list] = hr2hr_sparse(H_xyz,mode)
if nargin <2
    mode = 'num';
end
if strcmp(mode,'num')
    [NRPTS,~]=size(H_xyz);
    WAN_NUM=length(H_xyz(1).Hnum);
    Hnum_list{NRPTS} = sparse(WAN_NUM,WAN_NUM);
    %Hnum_list= zeros(WAN_NUM,WAN_NUM,NRPTS);
    vector_list = zeros(NRPTS,3);
   for i =1:NRPTS
        %Hnum_list(:,:,i) = H_xyz(i).Hnum./H_xyz(i).Degen;
        Hnum_list{i} = sparse(H_xyz(i).Hnum./double(H_xyz(i).Degen));
        vector_list(i,:) = H_xyz(i).vector;
   end
    hr_sparse.HnumL = Hnum_list;
    hr_sparse.vectorL = vector_list;
elseif strcmp(mode,'sym')
    
end
end