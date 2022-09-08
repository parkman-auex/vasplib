function H_hr = directsum_hr(H_hr1,H_hr2)
    [NRPTS,~]=size(H_hr1);
    H_hr = H_hr1;
    for i = 1:NRPTS
        H_hr(i).Hcoe = blkdiag(H_hr1(i).Hcoe,H_hr2(i).Hcoe);
        H_hr(i).Hnum = blkdiag(H_hr1(i).Hnum,H_hr2(i).Hnum);
    end
end