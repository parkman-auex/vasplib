function [dH_dkx_fun,dH_dky_fun,dH_dkz_fun] = Ham_diff(Ham_obj)
syms k_x k_y k_z real;
switch class(Ham_obj)
    case "Htrig"     
        symL = Ham_obj.HsymL_trig;
        dsym_dkx = diff(symL,k_x);
        dsym_dky = diff(symL,k_y);
        dsym_dkz = diff(symL,k_z);
    case "HK"
        symL = Ham_obj.HsymL;
        dsym_dkx = diff(symL,k_x);
        dsym_dky = diff(symL,k_y);
        dsym_dkz = diff(symL,k_z);
end
nbands = Ham_obj.Basis_num;
numL = reshape(Ham_obj.HnumL, nbands^2, []);
dH_dkx = reshape(numL*dsym_dkx.', nbands, nbands);
dH_dky = reshape(numL*dsym_dky.', nbands, nbands);
dH_dkz = reshape(numL*dsym_dkz.', nbands, nbands);

dH_dkx_fun = matlabFunction(dH_dkx,'Vars',[k_x,k_y,k_z]);
dH_dky_fun = matlabFunction(dH_dky,'Vars',[k_x,k_y,k_z]);
dH_dkz_fun = matlabFunction(dH_dkz,'Vars',[k_x,k_y,k_z]);
end