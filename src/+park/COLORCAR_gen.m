function [COLORCAR,WEIGHTCAR] = COLORCAR_gen(WAVECAR,HSVCAR)
    [~,norb,kn] = size(WAVECAR);
    temp.rgb = [0 0 0];
    COLORCAR = repmat(temp,norb,kn);
    WEIGHTCAR = zeros(norb,kn);
    for ki = 1:kn
        WAVECAR_one = WAVECAR(:,:,ki);
        for orbi = 1:norb
            WAVEFUNC = WAVECAR_one(:,orbi);
            rgb = COLOR_one_gen(WAVEFUNC,HSVCAR);
            COLORCAR(orbi,ki).rgb = rgb;
            hsv_temp = rgb2hsv(rgb);
            WEIGHTCAR(orbi,ki) = hsv_temp(1)*100+hsv_temp(3)*10+hsv_temp(2)*1;
        end
    end
end