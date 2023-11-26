function Ham_obj_B = add_zeeman_field(Ham_obj,B,opts)
% H_B = g*mu_B/hbar*s*B = mu_B*sigma*B (default: g=2, s=hbar/2*sigma)
arguments
    Ham_obj
    B double = [1 0 0] % Tesla, [Bx By Bz]
    opts.S double = 1/2
    opts.L double = 0
    opts.J double = 0
    opts.g_factor double = 2
    opts.spin_mode {mustBeMember(opts.spin_mode,{'new','old'})} = 'new'
    opts.natural_unit = false
end
%%
if isempty(Ham_obj.HnumL)
    error("Please input a numerical HR obj")
end
%%
sigma_x = [0 1;1 0]; sigma_y = [0 -1i;1i 0]; sigma_z = [1 0;0 -1];
if opts.natural_unit
    mu_B = 1;
else
    mu_B = 5.7884e-5; % eV/Tesla
end
s_dot_B = mu_B.*(B(1)*sigma_x + B(2)*sigma_y + B(3)*sigma_z);
%%
Nproj = Ham_obj.Nbands/2;
if opts.J == 0
    opts.J = opts.L + opts.S;
end
I_N_proj = diag(ones(1,Nproj) .* opts.J .* opts.g_factor);

if opts.spin_mode == "new"
    disp("spin indexes of projs = 1-up, 1-dn, 2-up, 2-dn, ...")
    disp("This is used in Wannier90 v2.x and newer")
    H_B = kron(I_N_proj, s_dot_B);
elseif opts.spin_mode == "old"
    disp("spin indexes of projs = 1-up, 2-up, ... 1-dn, 2-dn, ...")
    disp("This is used in Wannier90 v1.2 and Wannier Tools")
    H_B = kron(s_dot_B, I_N_proj);
end
Ham_obj_B = Ham_obj + H_B;
end