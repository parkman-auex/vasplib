%% Repeat  
% https://www.guanjihuan.com/archives/3932
% https://www.nature.com/articles/srep19018
% 10.1143/JPSJ.74.1674 
%% Chern 2 model
syms t1 t2 t3 k_x k_y m theta phi real;
%
% phi   = 0;
sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];
%%
d1 = Trig(2*t1*cos(k_x),sigma_x);
d2 = Trig(2*t1*cos(k_y),sigma_y);
d3 = Trig(m+2*t3*sin(k_x)+2*t3*sin(k_y)+2*t2*cos(k_x+k_y),sigma_z);
%
Type = 'sincos';
H_0 = Htrig(2,'Type',Type);
H_0 = H_0+d1+d2+d3;
%% para
t1 = 1;
t2 = 1;
t3 = 0.5;
m  = -1;
fermi = 0.5;
%% band
H_0n = H_0.Subsall();
H_0n = H_0n <'POSCAR_Chern2';
H_0n = H_0n <'KPOINTS_Chern2';
% [klist_l,kpoints_l,kpoints_name] = H_0n.kpath_information();
EIGENCAR = H_0n.EIGENCAR_gen();
bandplot(EIGENCAR,[-5,5],'title','Chern 2 bulk','Color','b','KPOINTS','KPOINTS_Chern2','POSCAR','POSCAR_Chern2');
%%
tic;
sizes = 51;
[BCCAR,Grid,~,klist] = H_0n.BC_2D( ...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-0.5,-0.5,0],'kdir1',[1,0,0],'kdir2',[0,1,0],...
    'plot',false,'oneshot',true);
[fig,ax] = vasplib.BCplot2D(BCCAR,Grid,H_0n.Rm,'BZ',true);
title('integral')
fprintf('Chern Number is %7.5f\n',sum(BCCAR,'all')/(2*pi));
toc;
%%
tic;
sizes = 51;
[klist_r,~,klist_r_plot,sizemesh,Gk_,Grid] = vasplib.kmesh2D(H_0n.Rm,...
    'knum1',sizes,'knum2',sizes, ...
    'kstart',[-0.5,-0.5,0],'kdir1',[1,0,0],'kdir2',[0,1,0]);
BC_test =  BC_kubo_formula(H_0n ,klist_r_plot);
BCCAR = reshape(BC_test(:,1),sizemesh);
% A3 = zeros(size(klist_r_plot));A3(:,3) = BC_test(:,1);
% [fig,ax] = BZplot(Haldane.Rm);
% vasplib_tool.quiverplot(klist_r,real(A3)*10,'b','Haldane:kubo num',fig,ax);
[fig,ax] = vasplib.BCplot2D(BCCAR,Grid,H_0n.Rm,'BZ',true);
toc;
title('kubo')
fprintf('Chern Number is %7.5f\n',sum(BCCAR,'all')/(2*pi));
% Hkfun =  H_0n.Hfun;
% Hkfun2D = @(k_x,k_y) Hkfun(k_x,k_y,0);
% %%
% knodes = 100;
% Rm = [1 0;0 1];
% [Bf,Bc_mat,fig] = BerryPhase_Discrete(Hkfun2D,Rm,knodes);



