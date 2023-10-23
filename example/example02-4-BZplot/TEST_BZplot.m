Rm = POSCAR_read('POSCAR');
Gk = eye(3)*2*pi/Rm;
Gk = Gk.';
%% bulk
[fig,ax] = vasplib.BZplot(Rm,'color','none');
%% Rm_xy = Rm(1:2,:)
SURFACE_CARD = [1 0 0;0 1 0;0 0 1];
Rm_xy = SURFACE_CARD* Rm;
Gk_tmp = (eye(3)*2*pi/Rm_xy).';
%Rm_xy(3,:) = [0 0 inf];
[fig,ax] = vasplib.BZplot(Rm_xy,'mode','2D','color','r','ax',ax,'label',false,...
'OriginPoint',Gk_tmp(3,:));
%% Rm_13 = Rm(1,3,:)
SURFACE_CARD = [1 0 0;-2/3 -1/3 1;-0.5 -1 0];
Rm_13 = SURFACE_CARD * Rm;
Gk_tmp = (eye(3)*2*pi/Rm_13).';
%Rm_xy(3,:) = [0 0 inf];
[fig,ax] = vasplib.BZplot(Rm_13,'mode','2D','color','g','ax',ax,'label',false,...
'OriginPoint',1*Gk_tmp(3,:));
%% Rm_23 = Rm(2,3,:)
SURFACE_CARD = [0.5 1 0;-2/3 -1/3 1;1 0 0];
Rm_23 =  SURFACE_CARD * Rm;
Gk_tmp = (eye(3)*2*pi/Rm_23).';
%Rm_xy(3,:) = [0 0 inf];
[fig,ax] = vasplib.BZplot(Rm_23,'mode','2D','color','b','ax',ax,'label',false,...
'OriginPoint',Gk_tmp(3,:));

Rm = POSCAR_read('POSCAR_166');
Gk = eye(3)*2*pi/Rm;
Gk = Gk.';
%% bulk
[fig,ax] = vasplib.BZplot(Rm,'color','none');
%% Rm_xy = Rm(1:2,:)
SURFACE_CARD = [1 0 0;0 1 0;0 0 1];
Rm_xy = SURFACE_CARD* Rm;
Gk_tmp = (eye(3)*2*pi/Rm_xy).';
%Rm_xy(3,:) = [0 0 inf];
[fig,ax] = vasplib.BZplot(Rm_xy,'mode','2D','color','r','ax',ax,'label',false,...
'OriginPoint',Gk_tmp(3,:));
%% Rm_13 = Rm(1,3,:)
SURFACE_CARD = [1 0 0;-2/3 -1/3 1;-0.5 -1 0];
Rm_13 = SURFACE_CARD * Rm;
Gk_tmp = (eye(3)*2*pi/Rm_13).';
%Rm_xy(3,:) = [0 0 inf];
[fig,ax] = vasplib.BZplot(Rm_13,'mode','2D','color','g','ax',ax,'label',false,...
'OriginPoint',1*Gk_tmp(3,:));
%% Rm_23 = Rm(2,3,:)
SURFACE_CARD = [0.5 1 0;-2/3 -1/3 1;1 0 0];
Rm_23 =  SURFACE_CARD * Rm;
Gk_tmp = (eye(3)*2*pi/Rm_23).';
%Rm_xy(3,:) = [0 0 inf];
[fig,ax] = vasplib.BZplot(Rm_23,'mode','2D','color','b','ax',ax,'label',false,...
'OriginPoint',Gk_tmp(3,:));