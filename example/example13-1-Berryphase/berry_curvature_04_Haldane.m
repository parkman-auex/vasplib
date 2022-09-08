% Haldane
syms h theta phi real;
%
% phi   = 0;
sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];
%%
syms a t t_2 M phi delta k_x k_y k_z real;
%
%
Rm = [1 0;-1/2 sqrt(3)/2]*a;
Haldane= Htrig(2);

t1 = [1/3,2/3]*a;
t2 = [2/3,1/3]*a;


e1 = t1-t2+[0  0]*Rm;
e2 = t1-t2+[1  0]*Rm;
e3 = t1-t2+[0 -1]*Rm;
d1 = [1 0]*Rm;
d2 = [0 1]*Rm;
d3 = [-1 -1]*Rm;

% H12 = -t*(...
%     cos(e1*[k_x;k_y]) +cos(e2*[k_x;k_y])+cos(e3*[k_x;k_y])...
%     - ...
%     1i*sin(e1*[k_x;k_y]) +sin(e2*[k_x;k_y])+cos(e3*[k_x;k_y]));
% 
% H21 = conj(H12);
H12trig_r = Trig(-t*(cos(e1*[k_x;k_y]) +cos(e2*[k_x;k_y])+cos(e3*[k_x;k_y])),sigma_x);
H12trig_i = Trig(-t*(sin(e1*[k_x;k_y]) +sin(e2*[k_x;k_y])+sin(e3*[k_x;k_y])),sigma_y);
H11_0 = Trig(-2*t_2*cos(phi)*(cos(d1*[k_x;k_y]) +cos(d2*[k_x;k_y])+cos(d3*[k_x;k_y])),sigma_0);
H11_z = Trig((M-2*t_2*sin(phi))*(sin(d1*[k_x;k_y]) +sin(d2*[k_x;k_y])+sin(d3*[k_x;k_y])),sigma_z);

a = 1;
%
Haldane= Haldane+ H12trig_r + H12trig_i+ H11_0+ H11_z;
Haldane.Rm = double(subs(Rm));


disp(Haldane);
%%
t   = double(1);
t_2  = 0.06;
M   =    0.0;
phi =    pi/2;
Haldane_n = Haldane.Subsall();
Haldane_n = Haldane_n < 'POSCAR';
Haldane_n = Haldane_n < 'KPOINTS';
%%
EIGENCAR = Haldane_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-3,3]);
%%
Haldane_hsym = subs(sym(Haldane));
%% BC from haldane definition
tic;
sizes = 51;
Bc = vasplib.BC_definition(Haldane_hsym,k_x,k_y);
Bcfun = matlabFunction(Bc,'Var',[k_x,k_y,k_z]);
[BCCAR,Grid,klist] = vasplib.BerryCuvature_fun(Bcfun, ...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-1,-1,0],'kdir1',[2,0,0],'kdir2',[0,2,0],...
    'plot',true,'Bcfun',true);
A = [zeros((sizes-1)*(sizes-1),2),real(BCCAR(:))];
[fig,ax] = BZplot(Haldane.Rm);
vasplib_tool.quiverplot(klist,real(A)*10,'r','Haldane:definition',fig,ax);
[fig,ax] = vasplib.BCplot2D(BCCAR,Grid,Haldane.Rm,'BZ',true);
toc;
%%  BC from haldane kubo
tic;
sizes = 51;
Bc = vasplib.BC_kubo_sym(Haldane_hsym,k_x,k_y);
Bcfun = matlabFunction(Bc,'Var',[k_x,k_y,k_z]);
[BCCAR,Grid,klist] = vasplib.BerryCuvature_fun(Bcfun, ...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-0.5,-0.5,0],'kdir1',[1,0,0],'kdir2',[0,1,0],...
    'plot',true,'Bcfun',true);
A = [zeros((sizes-1)*(sizes-1),2),real(BCCAR(:))];
[fig,ax] = BZplot(Haldane.Rm);
vasplib_tool.quiverplot(klist,real(A)*10,'r','Haldane:kubo',fig,ax);
[fig,ax] = vasplib.BCplot2D(BCCAR,Grid,Haldane.Rm,'BZ',true);
toc;
%% BC from haldane kubo semi-numerical Yang Fan
tic;
sizes = 51;
[klist_r,~,klist_r_plot,sizemesh,Gk_,Grid] = vasplib.kmesh2D(Haldane.Rm,...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-1,-1,0],'kdir1',[2,0,0],'kdir2',[0,2,0]);
BC_test =  BC_kubo_formula(Haldane_n,klist_r_plot);
BCCAR = reshape(BC_test(:,1),sizemesh);
A3 = zeros(size(klist_r_plot));A3(:,3) = BC_test(:,1);
[fig,ax] = BZplot(Haldane.Rm);
vasplib_tool.quiverplot(klist_r,real(A3)*10,'b','Haldane:kubo num',fig,ax);
[fig,ax] = vasplib.BCplot2D(BCCAR,Grid,Haldane.Rm,'BZ',true);
toc;
%% BC from haldane kubo semi-numerical Yang Fan 

%% To be tested
test = 1; 
if test == 2
%% BF from haldane definition
[kpoints,nodes,kpoints_name] = KPOINTS_read('KPOINTS_integral');
[kloop,~,~,~,~]= kpathgen3D(Haldane.Rm,kpoints,nodes,kpoints_name);
iBF_sym  = vasplib.BerryPhaseLine_definition_sym_i(Haldane_hsym);
iBF_fun = matlabFunction(iBF_sym,'Var',[sym('k_x'),sym('k_y'),sym('k_z'),sym('delta_k_x'),sym('delta_k_y'),sym('delta_k_z')]);
BF1 = vasplib.BerryPhaseLine_fun(iBF_fun,kloop,'Rm',Haldane.Rm,'plot',true,'funtype','Ffun');
%% BD from haldane parallel-transport gauge discrete symbolic
[kpoints,nodes,kpoints_name] = KPOINTS_read('KPOINTS_integral');
[kloop,~,~,~,~]= kpathgen3D(Haldane.Rm,kpoints,nodes,kpoints_name);
[W,~] =eig(Haldane_hsym);
WAVE_fun = matlabFunction(vasplib.NomalizeEigenvector(W(:,1)),'Var',[sym('k_x'),sym('k_y'),sym('k_z')]);
BF2 = vasplib.BerryPhaseLine_fun(WAVE_fun,kloop,'Rm',Haldane.Rm,'plot',true,'funtype','Wfun');
end