%% init_hr
function H_hr = init_hr(wan_num)
    H_hr_ =struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
    H_hr(1)  = H_hr_;
    H_hr(1).seq = 1;
    H_hr(1).Degen = 1;
    H_hr(1).vector = [0,0,0];
    H_hr(1).Hcoe = sym(zeros(wan_num));
    H_hr(1).Hnum = zeros(wan_num);
    %H_hr(1).Hsym = sym(zeros(wan_num));
end