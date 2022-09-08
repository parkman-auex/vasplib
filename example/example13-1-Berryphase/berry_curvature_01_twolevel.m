%% berry_curvature_01_twolevel
syms h theta phi real;
%
% phi   = 0;
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];
%
h_vector(1) = h*sin(theta)*cos(phi);
h_vector(2) = h*sin(theta)*sin(phi);
h_vector(3) = h*cos(theta);

H = sigma_x*h_vector(1)+sigma_y*h_vector(2)+sigma_z*h_vector(3);
% first WaveFunc
[WAVECAR,EIGENCAR_diag] = eig(H);
EIGENCAR = diag(simplify(EIGENCAR_diag));
u1 = simplify(vasplib.NomalizeEigenvector(WAVECAR(:,2))); 
% the second is -h eigenstate
u2 = simplify(vasplib.NomalizeEigenvector(WAVECAR(:,1))); 
WAVECAR = [u1,u2];
% in 10.1103/RevModPhys.82.1959
u_minus = [sin(theta/2)*exp(-1i*phi);-cos(theta/2)];
u_plus = [cos(theta/2)*exp(-1i*phi);sin(theta/2)];
WAVECAR2 = [u_plus,u_minus];
%%
fprintf('Use eigenstates by matlab:\n');
disp(simplify(WAVECAR\H*WAVECAR));
disp(WAVECAR);               % latex()
fprintf('Use eigenstates by Article:\n');
disp(simplify(WAVECAR2\H*WAVECAR2));% latex()
disp(WAVECAR2);
%%
Berry_curvature = vasplib.Berry_curvature_D2(u1,theta,phi);
 fprintf('Berry_curvature : %s !\n',string(Berry_curvature));
%%
fprintf('add a random phase factor.\n');
syms ksi real;
u1_prime = u1*exp(1i*ksi);
Berry_curvature_gaugefree = vasplib.Berry_curvature_D2(u1_prime,theta,phi);
 fprintf('Gauge free: Berry_curvature '' : %s !\n',...
     string(Berry_curvature_gaugefree));
%%
% syms theta phi h x y z real;
% h = (x^2+y^2+z^2)^(1/2);
% theta = acos(z/h);
% phi = acos(x/(h*sin(theta)));
% J = simplify(jacobian([phi,cos(theta)],[y,z]));
% disp(latex(0.5*det(J)));

