%% useful_tools
sigma_0 = pauli_matric(0); sigma_x = pauli_matric(1); sigma_y = pauli_matric(2); sigma_z = pauli_matric(3);
tau_0   = pauli_matric(0); tau_x   = pauli_matric(1); tau_y   = pauli_matric(2); tau_z   = pauli_matric(3);
%% 
syms Delta m v lambda u1 u2;
syms k_x k_y k_z d real ;
H_2D = HK(4,2);
%% 
% % the coeffs is not allowed use x y z

H_2D = H_2D ...
    + Term( -Delta              ,    sigma_z*tau_0)...
    + Term( v*k_x,                   sigma_x*tau_x)...
    + Term( v*k_y,                   sigma_x*tau_y);

%     + Term( lambda*(k_x^2-k_y^2),    sigma_z*tau_x)...
%     + Term( 2*lambda*(k_x * k_y),    sigma_z*tau_y)...
    
H_3D = H_2D ...
    + Trig( u1*(cos(d*k_z) )         ,    sigma_z*tau_0)...
    + Trig( -u2*(sin(d*k_z))          ,    sigma_y*tau_z)...
    ;
H_3D = H_3D < 'POSCAR';
H_3D = H_3D < 'KPOINTS';
%% 
% % para

b = 1.424e-10;
h_bar =6.582119514e-16; %eV⋅s
v = 7e5 * h_bar /b  ; % m/s
%% 
d =3.3;
% 

%v =0;
u1 =0.233;
u2 =0.23;

Delta =0.213;
tolerance = 0.01;
kmesh = [10,10,0.825525678548721];
H_3Dn = H_3D.Subsall();
%
EIGENCAR = H_3Dn.EIGENCAR_gen();
bandplot(EIGENCAR,[-1,1]);

%[k_list_cart,klist_f,~]=H_3Dn.findnodes(kmesh,2,tolerance);
H_3Dn = H_3Dn < 'KPOINTS_findnode';
import vasplib_tool.*;
[fig,ax] = creat_figure();
kzf = linspace(0,1,100);
Gk = (eye(3)*2*pi)/H_3Dn.Rm;
for i = 1:length(kzf )
    H_3Dn.klist_r(:,3) = kzf(i)*(Gk(3,3));
    EIGENCAR = H_3Dn.EIGENCAR_gen();
    [bandgap,label] = min(EIGENCAR(3,:)- EIGENCAR(2,:));
    scatter(i,EIGENCAR(2,label));
    scatter(i,EIGENCAR(3,label));
    if bandgap < tolerance
        fprintf('find it: %7.4f eV (%6.3f, %6.3f,%6.3f)\n',bandgap,H_3Dn.klist_s(i,1),H_3Dn.klist_s(i,2),kzf(i));
    else
        fprintf('min gap: %7.4f eV (%6.3f, %6.3f,%6.3f)\n',bandgap,H_3Dn.klist_s(i,1),H_3Dn.klist_s(i,2),kzf(i));
    end

end