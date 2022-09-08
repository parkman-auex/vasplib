function OUT_H_xyz = zero_H_xyz(H_xyz,WAN_NUM)
    OUT_H_xyz = H_xyz;
    [NRPTS,~]=size(H_xyz);
    for i =1:NRPTS
        % not good
        OUT_H_xyz(i).Hstr = zeros(WAN_NUM);
        OUT_H_xyz(i).Hsym = zeros(WAN_NUM);
        OUT_H_xyz(i).Hcoe = zeros(WAN_NUM);
        % not good
        OUT_H_xyz(i).Hnum = zeros(WAN_NUM);
    end
end