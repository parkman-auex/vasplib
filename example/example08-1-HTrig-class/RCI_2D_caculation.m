%% RCI_2D
% RCI model
%% RCI model
% Useful matrix

clear;
Gamma_1_4 = sym(gamma_matric(1,4));
Gamma_1_5 = sym(gamma_matric(1,5));
Gamma_2_4 = sym(gamma_matric(2,4));
Gamma_2_5 = sym(gamma_matric(2,5));
%% useful tool


% k·p model

syms M u real;
syms k_x k_y k_z real;
assume(k_x <= pi & k_x >= -pi);
assume(k_y <= pi & k_y >= -pi);
assume(k_z <= pi & k_z >= -pi);
U = expm(-1i*pi*sym(gamma_matric(2,5))/4);
%
% M       = M0-M2*(k_x^2+k_y^2)-M1*(k_z^2);
% E0k     = C0+C2*(k_x^2+k_y^2)+C1*(k_z^2);
% k_plus  = k_x + 1i* k_y;
% k_minus = k_x - 1i* k_y;
%
%%
RCI = Htrig(4)
%
RCI = RCI ...
    +Trig(sin(k_x) ,  gamma_matric(1) )...
    +Trig(sin(k_y) ,  gamma_matric(2))...
    +Trig(M-cos(k_x)-cos(k_y)   , gamma_matric(3) )...
    ;
RCI_sym = sym(RCI);
syms m_1 m_2 m_3 m_4;
syms m_0 theta_0 theta_x theta_y m_x m_y  real;
m_1 = m_y*sin(theta_y);
m_3 = m_x*sin(theta_x);
m_2 = m_y*cos(theta_y);
m_4 = m_x*cos(theta_x);
m_x = m_0*cos(theta_0);
m_y = m_0*sin(theta_0) ;
delta_H = Htrig(4,Trig(m_1,Gamma_1_4)+Trig(m_2,Gamma_2_4)+Trig(m_3,Gamma_1_5)+Trig(m_4,Gamma_2_5));
delta_H.Htrig_sym   
RCI_PT = RCI + delta_H;
RCI_PT = RCI_PT <= 'POSCAR';
RCI_PT = RCI_PT <'KPOINTS';
%%
syms theta_x theta_y theta_0 m_0 real
RCIPT = RCI_PT.Htrig_sym
C =m_1*gamma_matric(4)+m_3*gamma_matric(5)
m_x = m_0*cos(theta_0);
m_y = m_0*sin(theta_0) ;
% theta_x = theta_y;
simplify(subs(C*RCIPT+RCIPT*C))
% Bulk band
% para
% band

M = 1;
m_0 = 0.5;
theta_0 = 0.5*pi;
theta_x = 0.2*pi;
theta_y = 0.3*pi;
calculate =1;
if calculate ==1
RCI_n = RCI_PT.Subsall();
% bulk
EIGENCAR_disk = RCI_n.EIGENCAR_gen();
[klist_l,kpoints_l,kpoints_name] = RCI_n.kpath_information();
bandplot(EIGENCAR_disk,[-3,3],klist_l,kpoints_l,kpoints_name, ...
    'Color','r','title','');
title({['RCI_{PT}','M=1,m=0.5'];['\theta_0 = ',char(string(theta_0/pi))...
,'\pi'];[ '\theta_x = ' ,char(string(theta_x/pi))...
,'\pi'];[' \theta_y = ' ,char(string(theta_y/pi)),'\pi']})
end
%%
calculate =2;
if calculate ==2
% clear('theta_x');clear('theta_y')[sym('theta_x'),sym('theta_y')]
clear('m_0');
RCI_PT.HcoeL = subs(RCI_PT.HcoeL);
RCI_n = RCI_PT.Subsall('para',sym('m_0'));
EIGENCAR_cell = RCI_n.EIGENCAR_gen('paraname','m_0','para',(0:0.02:2).');
bandplot(EIGENCAR_cell,[-3,3],klist_l,kpoints_l,kpoints_name,'title','RCIPT-m_0/M = 0:2 ','Color',@parula);
colorbar(gca,'Ticks',[0 0.25 0.5 0.75 1],'TickLabels',{'0','1/2','1','3/2','2'});
end

%% slab

calculate =3;
M =1.5;
m_0 = 0;
theta_0 = 0*pi;
syms theta_x theta_y real;
if calculate ==3
%RCI_d = RCI_PT;
RCI_d = RCI.discretize([30 0 0]);
RCI_d.Htrig_sym
RCI_d = RCI_d <'KPOINTS_slaby';
[klist_l,kpoints_l,kpoints_name] = RCI_d.kpath_information();
RCI_d.HcoeL = subs(RCI_d.HcoeL);
% RCI_d = RCI_d.Subsall('para',[sym('theta_x'),sym('theta_y')]);
% [EIGENCAR_slab] = RCI_d.EIGENCAR_gen_slab('norb',-1,'paraname',[sym('theta_x'),sym('theta_y')],'para',[(0:0.1:1).',-(0:0.1:1).']);
theta_x = pi/4;
theta_y = pi/4;
RCI_d = RCI_d.Subsall('slab');
[EIGENCAR_slab] = RCI_d.EIGENCAR_gen_slab();
bandplot(EIGENCAR_slab,[-2,2],klist_l,kpoints_l,kpoints_name,'title','RCI-slab','Color',@jet);
end
%bandplot(EIGENCAR_slab,[-2,2],klist_l,kpoints_l,kpoints_name,'title','RCI-slab','Color',@jet);
% LDOC in CORNER
% 
% 
% 
% Corner
% bugs here

%%
calculate =4;
M =1;
m_0 = 0.2;
theta_0 = 0*pi;
RCI_d = RCI;
RCI_d = RCI_d.discretize([20 20 0])
RCI_d = RCI_d <'KPOINTS_wire';
RCI_d = RCI_d.Subsall('disk');
[EIGENCAR_disk,...
WAVECAR_disk] = RCI_d.EIGENCAR_gen_wire('klist',[0,0,0]);
figure();
plot(EIGENCAR_disk,'-o')
PARCHG_gen(RCI_d.orbL,WAVECAR_disk(:,800:801)*3)

