%% for
function MAGMOM=MAGMOM_gen(mag_list,MagAmp,mode,ions_num)
if nargin <3
    mode = 'from_maglist';
end
mags_num = length(mag_list);
if nargin <4
    ions_num = mags_num;
end
MAGMOM = "";
switch mode
    case 'from_maglist'
        
    case 'PM'
        for i = 1:mags_num
            MAGMOM = MAGMOM + string("0"+" ");
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num)+"*"+0;
        end
    case 'FM'
        for i = 1:mags_num
            MAGMOM = MAGMOM + string(sign(mag_list(i))*MagAmp)+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num)+"*"+0;
        end
    case 'FMx'
        for i = 1:mags_num
            MAGMOM = MAGMOM + string(MagAmp)+" 0 0"+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num*3)+"*"+0;
        end
    case 'FMy'
        for i = 1:mags_num
            MAGMOM = MAGMOM + "0 "+string(MagAmp)+" 0"+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num*3)+"*"+0;
        end
    case 'FMz'
        for i = 1:mags_num
            MAGMOM = MAGMOM +"0 0 "+string(MagAmp)+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num*3)+"*"+0;
        end
    case 'AFM'
        for i = 1:mags_num
            MAGMOM = MAGMOM + string(sign(mag_list(i))*MagAmp)+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num)+"*"+0;
        end
    case 'AFMx'
        for i = 1:mags_num
            MAGMOM = MAGMOM + string(sign(mag_list(i))*MagAmp)+" 0 0"+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num*3)+"*"+0;
        end
    case 'AFMy'
        for i = 1:mags_num
            MAGMOM = MAGMOM + "0 "+string(sign(mag_list(i))*MagAmp)+" 0"+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num*3 )+"*"+0;
        end
    case 'AFMz'
        for i = 1:mags_num
            MAGMOM = MAGMOM +"0 0 "+string(sign(mag_list(i))*MagAmp)+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num*3)+"*"+0;
        end
    case '120'
        disp('suppose a^b = \pi/6');
        for i = 1:mags_num
            switch mag_list(i)
                case 1
                    kep_string = string(MagAmp)+" 0 0";
                case 2
                    kep_string = "0 "+string(-MagAmp)+" 0";
                case 3
                    kep_string = string(-MagAmp)+" "+string(MagAmp)+" 0";
            end
            MAGMOM = MAGMOM +kep_string+" ";
        end
        nomag_num = ions_num-mags_num;
        if nomag_num >0
            MAGMOM = MAGMOM + string(nomag_num*3)+"*"+0;
        end        
end
end
