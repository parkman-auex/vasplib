%% hr_sparse2hr
% we will translate a 3D mat a lisr_for NRPT list(V) and a list for Degen()
% suppose all Degen = 1
function [H_xyz] = hr_sparse2hr(hr_sparse)
    Hnum_list = hr_sparse.HnumL ;
    vector_list = hr_sparse.vectorL ;
    [NRPTS,~]=size(vector_list);
    WAN_NUM=length(Hnum_list(:,:,1));
    
    % H_xyz[]
    H_xyz_ =struct('seq',[],'vector',[0, 0 ,0],'Degen',[],'Hcoe',zeros(WAN_NUM),'Hnum',zeros(WAN_NUM));
    H_xyz = repmat(H_xyz_ ,[NRPTS,1]);    % Hamiltonian of every cell;
    
    for i =1:NRPTS
        H_xyz(i).seq = i ;
        H_xyz(i).Hnum = full(Hnum_list{i});
        H_xyz(i).Hcoe = full(Hnum_list{i}) ;
        H_xyz(i).Degen = 1;
        H_xyz(i).vector = vector_list(i,:);
    end
end