function H_E = add_electric_field(Ham_obj,E)
arguments
    Ham_obj
    E double = [1 0 0] % V/Ang, [Ex Ey Ez]
end
orbL_cart = Ham_obj.orbL * Ham_obj.Rm;
E_onsite = orbL_cart * E';
H_E = diag(E_onsite);
end