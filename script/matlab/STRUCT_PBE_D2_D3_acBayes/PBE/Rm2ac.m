function [a,c]  =Rm2ac(Rm,lattice)
    switch lattice
        case 'hex'
            a = Rm(1,1);
            c = Rm(3,3);
        case 'tetra'
            a = Rm(1,1);
            c = Rm(3,3);
        case ''
    end
end