function COLOR_one = COLOR_one_gen(WF,HSVCAR)
    n = abs(WF.*conj(WF));
    %[WAN_NUM,~] = size(HSVCAR);
    %HSVCAR(:,2:3) = ones(WAN_NUM,2);
    RGBCAR =  hsv2rgb(HSVCAR);
    COLOR_one = sum(RGBCAR.*n,1);
    COLOR_one = normalize(COLOR_one,'range');
end