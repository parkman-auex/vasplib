%% reshape_hr
function H_xyz = reshape_hr(H_xyz,trans_orb)
[NRPTS,~]=size(H_xyz);
WAN_NUM=length(H_xyz(1).Hnum);
trans_orb1 =  trans_orb(1);
trans_orb2 =  trans_orb(2);

trans_matrix = eye(WAN_NUM);
trans_matrix(trans_orb1 ,trans_orb2) = 1 ; 
trans_matrix(trans_orb2 ,trans_orb1) = 1 ; 
trans_matrix(trans_orb1 ,trans_orb1) = 0 ; 
trans_matrix(trans_orb2 ,trans_orb2) = 0 ; 
for i = 1:NRPTS
    H_xyz(i).Hcoe = trans_matrix * H_xyz(i).Hcoe * trans_matrix' ;
    H_xyz(i).Hnum = trans_matrix * H_xyz(i).Hnum * trans_matrix' ;
end

end