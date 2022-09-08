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
a = 1;
t   = double(1);
t_2  = 0.06;
M   =    0.0;
phi =    sym(pi/2);
Rm = sym([1 0;-1/2 sqrt(3)/2]);
Type = 'list';
Haldane= Htrig(2,'Type',Type);

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


%
Haldane= Haldane+ H12trig_r + H12trig_i+ H11_0+ H11_z;
Haldane.Rm = double(subs(Rm));
%Haldane = Haldane.rewrite();

disp(Haldane);
%%

Haldane_n = Haldane.Subsall();
Haldane_n = Haldane_n < 'POSCAR';
Haldane_n = Haldane_n < 'KPOINTS';
%%
EIGENCAR = Haldane_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-3,3]);
%%
Haldane_hsym = subs(sym(Haldane));
%%
sizes = 51;
tic;
[BCCAR,Grid,~,klist] = Haldane_n.BC_2D( ...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-1,-1,0],'kdir1',[2,0,0],'kdir2',[0,2,0],...
    'plot',false,'oneshot',true,'BAND_index',1);
[fig,ax] = vasplib.BCplot2D(BCCAR,Grid,Haldane_n.Rm,'BZ',true,'shading',true);
fprintf('Chern Number is %7.5f\n',sum(BCCAR,'all')/(2*pi));
% %%  BC from haldane kubo
% tic;
% sizes = 51;
% Bc = vasplib.BC_kubo_sym(Haldane_hsym,k_x,k_y);
% Bcfun = matlabFunction(Bc,'Var',[k_x,k_y,k_z]);
% [BCCAR,Grid,klist] = vasplib.BerryCuvature_fun(Bcfun, ...
%     'knum1',sizes,'knum2',sizes, ...
%     'kstart',[-1,-1,0],'kdir1',[2,0,0],'kdir2',[0,2,0],...
%     'plot',true,'Bcfun',true);
% [fig,ax] = vasplib.BCplot2D(BCCAR,Grid,Haldane.Rm,'BZ',true);
toc;
%% BC from haldane kubo semi-numerical Yang Fan test for list Htrig
tic;
sizes = 51;
[klist_r,~,klist_r_plot,sizemesh,Gk_,Grid] = vasplib.kmesh2D(Haldane.Rm,...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-1,-1,0],'kdir1',[2,0,0],'kdir2',[0,2,0]);
BC_test =  BC_kubo_formula(Haldane_n,klist_r_plot);
BCCAR = reshape(BC_test(:,1),sizemesh);
% A3 = zeros(size(klist_r_plot));A3(:,3) = BC_test(:,1);
% [fig,ax] = BZplot(Haldane.Rm);
% vasplib_tool.quiverplot(klist_r,real(A3)*10,'b','Haldane:kubo num',fig,ax);
[fig,ax] = vasplib.BCplot2D(BCCAR,Grid,Haldane.Rm,'BZ',true,'shading',true);
toc;