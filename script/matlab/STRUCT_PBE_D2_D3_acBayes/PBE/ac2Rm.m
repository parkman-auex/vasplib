function Rm  =ac2Rm(a,c,lattice)
    switch lattice
        case 'hex'
            Rm = [a,0,0;1/2*a sqrt(3)*a/2 0;0 0 c];
        case 'tetra'
            Rm = [a,0,0;1/2*a sqrt(3)*a/2 0;0 0 c];
        case ''
    end
end