function orb_sym = Ymlsym(l,m,orb_symbk)
    if l == 1i && m == 1i
        orb_sym = orb_symbk;
    end
    switch l
        case 0
            switch m
                case 0
                    orb_sym = str2sym('1');
                otherwise
                    warning('check your input POSCAR');
            end
        case 1
            switch m
                case 0
                    orb_sym = str2sym('z');
                case -1
                    orb_sym = str2sym('y');
                case 1
                    orb_sym = str2sym('x');
                otherwise
                    warning('check your input POSCAR');
            end
        case 2
            switch m
                case 0
                    orb_sym = str2sym('z^2');
                case -1
                    orb_sym = str2sym('y*z');
                case 1
                    orb_sym = str2sym('x*z');
                case -2
                    orb_sym = str2sym('x*y');
                case 2
                    orb_sym = str2sym('x^2-y^2');
                otherwise
                    warning('check your input POSCAR');
            end
        case 3
            switch m
                case 0
                    orb_sym = str2sym('z^3');
                case -1
                    orb_sym = str2sym('y*z^2');
                case 1
                    orb_sym = str2sym('x*z^2');
                case -2
                    orb_sym = str2sym('x*y*z');
                case 2
                    orb_sym = str2sym('(x^2-y^2)*z');
                case -3
                    orb_sym = str2sym('y*(3*x^2-y^2)');
                case 3
                    orb_sym = str2sym('x*(x^2-3*y^2)');
                otherwise
                    warning('check your input POSCAR');
            end
        case -1
            switch m
                case 1
                    orb_sym = str2sym('1+x');
                case 2
                    orb_sym = str2sym('1-x');
                otherwise
                    warning('check your input POSCAR');
            end 
        case -2
            switch m
                case 1
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x+2^(-1/2)*y');
                case 2
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x-2^(-1/2)*y');
                case 3
                    orb_sym = str2sym('3^(-1/2)+6^(-1/2)*x+6^(-1/2)*x');                 
                otherwise
                    warning('check your input POSCAR');
            end
        case -3
            switch m
                case 1
                    orb_sym = str2sym('1+x+y+z');
                case 2
                    orb_sym = str2sym('1+x-y-z');
                case 3
                    orb_sym = str2sym('1-x+y-z'); 
                case 4
                    orb_sym = str2sym('1-x-y+z');
                otherwise
                    warning('check your input POSCAR');
            end
        case -4 % 
            switch m
                case 1
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x+2^(-1/2)*y');
                case 2
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x-2^(-1/2)*y');
                case 3
                    orb_sym = str2sym('3^(-1/2)+6^(-1/2)*x+6^(-1/2)*x');
                case 4
                    orb_sym = str2sym('2^(-1/2)*z+2^(-1/2)*z^2');
                case 5
                    orb_sym = str2sym('-2^(-1/2)*z+2^(-1/2)*z^2');
                otherwise
                    warning('check your input POSCAR');
            end
        case -5
            switch m
                case 1
                    orb_sym = str2sym('6^(-1/2)-2^(-1/2)*x-12^(-1/2)*z^2+2^(-1)*(x^2-y^2)');
                case 2
                    orb_sym = str2sym('6^(-1/2)+2^(-1/2)*x-12^(-1/2)*z^2+2^(-1)*(x^2-y^2)');
                case 3
                    orb_sym = str2sym('6^(-1/2)-2^(-1/2)*x-12^(-1/2)*z^2-2^(-1)*(x^2-y^2)');
                case 4
                    orb_sym = str2sym('6^(-1/2)+2^(-1/2)*x-12^(-1/2)*z^2-2^(-1)*(x^2-y^2)');
                case 5
                    orb_sym = str2sym('6^(-1/2)-2^(-1/2)*z+3^(-1)*(z^2)');
                case 6
                    orb_sym = str2sym('6^(-1/2)+2^(-1/2)*z+3^(-1)*(z^2)');
                otherwise
                    warning('check your input POSCAR');
            end
        otherwise
            warning('check your input POSCAR');
    end
end