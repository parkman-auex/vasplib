function add_Peierls_substitution(Ham_obj)
% H_B = g*mu_B/hbar*s*B = mu_B*sigma*B (default: g=2, s=hbar/2*sigma)
arguments
    Ham_obj HR
end
%%
if ~Ham_obj.num
    error("Please input a numerical HR obj")
end
%%
vectorL_cart = Ham_obj.vectorL * Ham_obj.Rm;
orbL_cart = Ham_obj.orbL * Ham_obj.Rm;
A = [1 0 0]; % B=1T
phi0 = 4.13567e-15; % Wb, T*m^2
phi0 = phi0 * 1e20; % Wb, T*Ang^2
%%
for i = 1:Ham_obj.WAN_NUM
    for j = 1:Ham_obj.WAN_NUM
        for k = 1:Ham_obj.NRPTS
            xyz_mid = ( orbL_cart(i,:) + orbL_cart(j,:) + vectorL_cart(k,:) )./2; % Ang
            
            A_field = A * xyz_mid'; % (y*B, 0 , 0)
            A_dl = A_field * ( orbL_cart(j,1) + vectorL_cart(k,1) - orbL_cart(i,1) ); % x_j - x_i

            Ham_obj.HnumL(i,j,k) = Ham_obj.HnumL(i,j,k) * exp(1i * 2*pi/phi0 * A_dl);
        end
    end
end

end