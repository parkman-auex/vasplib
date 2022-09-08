function Carbon_Nanotube_gen(m,n,mcell,ncell,opts)
arguments
    m double = 1;
    n double = 0;
    mcell double = 20;
    ncell double = 1;
    opts.Tube_Type {mustBeMember(opts.Tube_Type,{'Normal','Mobius'})} = 'Normal';
    opts.Lat_Type {mustBeMember(opts.Lat_Type,{'Square','Honeycomb'})} = 'Honeycomb';
    opts.POSCAR string = "POSCAR"
end
disp('遵循Chiral碳纳米管标准范式，请采用基矢夹角为锐角的石墨烯POSCAR');
disp('X为径向，Y为轴向，请通过调节手性数m、n来调节边界');
disp('Mobius模式下，建议ncell = 1，以限制碳纳米管宽度，减轻内圈结构畸变')
[Rm,sites,Atom_name,Atom_num,~,~] = POSCAR_readin(opts.POSCAR);
%% primitive cell of graphene
% Rm = [2.4560000896         0.0000000000         0.0000000000;
%       1.2280000448         2.1269584693         0.0000000000;
%       0.0000000000         0.0000000000        20.0000000000];
% Atom_name = "C";
% Atom_num = 2;
% sites(1) = struct('rc1',1/3,'rc2',1/3,'rc3', 1/4, 'name',Atom_name);
% sites(2) = struct('rc1',2/3,'rc2',2/3,'rc3', 1/4, 'name',Atom_name);
%% to rectangular ribbon
switch opts.Lat_Type
    case 'Honeycomb'
        L = gcd(m,n);
        Rotate1 = [m/L n/L 0; -(m+2*n)/L (2*m+n)/L 0; 0 0 1];
    case 'Square'
        Rotate1 = [1 0 0; 0 1 0; 0 0 1];
end
[Rm,sites] = supercell(Rotate1,Rm,sites,Atom_name,Atom_num,[0 0 0],'POSCAR_ribbon');
Atom_num = size(sites,2);
%% enforce diagonal Rm
Rm = diag([sqrt(Rm(1,1)^2+Rm(1,2)^2), sqrt(Rm(2,1)^2+Rm(2,2)^2), Rm(3,3)]);

Rotate2 = [mcell 0 0; 0 ncell 0; 0 0 1];
[Rm_ribbon,sites_ribbon] = supercell(Rotate2,Rm,sites,Atom_name,Atom_num,[0 0 0],'POSCAR_ribbon');
%% to nanotube
radius = Rm_ribbon(1,1)/2/pi;
box_radius = Rm_ribbon(1,1)/2/pi + Rm_ribbon(2,2)/2 + 2;
Rm_tube = diag(ones(1,3) .* box_radius * 2);

Atom_num = size(sites_ribbon,2);
for i = 1:Atom_num
    xyz_cart = [sites_ribbon(i).rc1 sites_ribbon(i).rc2 0]*Rm_ribbon;
    
    theta = xyz_cart(1)/radius;
    if strcmp(opts.Tube_Type,'Normal')
        phi = 0;
    elseif strcmp(opts.Tube_Type,'Mobius')
        phi = theta/2;
    end
    
    wide = xyz_cart(2) - Rm_ribbon(2,2)/2;
    z_cart = wide * cos(phi);
    radius_new = radius + wide * sin(phi);
    
    x_cart = cos(theta) *  radius_new + box_radius;
    y_cart = sin(theta) *  radius_new + box_radius;
    z_cart = z_cart + box_radius;
    
    xyz_frac = [x_cart, y_cart, z_cart] / Rm_tube;
    sites_tube(i) = struct('rc1',xyz_frac(1),'rc2',xyz_frac(2),'rc3',xyz_frac(3),'name',Atom_name);
end
POSCAR_gen(Rm_tube,sites_tube,Atom_name,Atom_num,'POSCAR_tube',-8);