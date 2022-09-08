% Haldane
syms h theta phi real;
%
% phi   = 0;
sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];
%%
syms a t t_2 M phi delta k_x k_y real;
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
t_2  = 0.2;
M   =    0.1;
%phi =    pi/2;
phi = pi/2;
Haldane_n = Haldane.Subsall();
Haldane_n = Haldane_n < 'POSCAR';
Haldane_n = Haldane_n < 'KPOINTS';
%%
EIGENCAR = Haldane_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-3,3]);
%%
sizes = 51;
Fig = Figs(1,2);
[BCCAR,Grid,~,klist] = Haldane_n.BC_2D( ...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-0.5,-0.5,0],'kdir1',[1,0,0],'kdir2',[0,1,0],...
    'plot',false,'oneshot',true);
vasplib_plot.BCplot2D(BCCAR,klist,Haldane_n.Rm,'BZ',true,'BZmode','2D','ax',Fig.axes(1));title(Fig.axes(1),'Calculation');
[klist_expand,BCCAR_expand] = kshift(klist,[-1,-0.5,0;2,-1,0;0,2,0],BCCAR,'cart',true,'Rm',Haldane_n.Rm);

vasplib_plot.BCplot2D(BCCAR_expand,klist_expand,Haldane_n.Rm,'BZ',true,'BZmode','2D','ax',Fig.axes(2));title(Fig.axes(2),'AutoExpand');
fprintf('Chern Number is %7.5f\n',sum(BCCAR,'all')/(2*pi));
%%
[fig,ax] = BZplot(Haldane.Rm);
A = [zeros((sizes-1)*(sizes-1),2),real(BCCAR(:))];
vasplib_plot.quiverplot(klist,real(A)*10,'r','Haldane:integral',ax);
