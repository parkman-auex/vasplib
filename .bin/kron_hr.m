%% kron_hr
function H_xyz = kron_hr(matrix,H_xyz)


[NRPTS,~]=size(H_xyz);
% WAN_NUM=length(H_xyz(1).Hnum);

for i = 1:NRPTS
    H_xyz(i).Hcoe = kron(matrix,H_xyz(i).Hcoe ) ;
    H_xyz(i).Hnum = kron(matrix,H_xyz(i).Hnum ) ;
end

end