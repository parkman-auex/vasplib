%% HRclass Tutorial
% 基于HRclass实现石墨烯model      
% 
% 2021.08
%% 
% * Author: parkman
% * Email：parkman@buaa.edu.cn
% Prepare

syms t real;
Grpahene_TB = HR(2);
A_vectorL = [0,0,0;1,0,0;0,-1,0];
B_vectorL = [0,0,0;-1,0,0;0,1,0];
%%
Grpahene_TB = Grpahene_TB.set_hop(t,1,2,A_vectorL,'sym');
Grpahene_TB = Grpahene_TB.set_hop(t,2,1,B_vectorL,'sym')
%%
Grpahene_TB.printout;
%%
Grpahene_TB_list = Grpahene_TB.rewrite()
%%
Grpahene_TB_list.list()
%%
Rm = [1,0,0;-0.5,sqrt(3)/2,0;0,0,1];
Grpahene_TB_list = Grpahene_TB_list.input_Rm(Rm);
orbL = [2/3,1/3,0;1/3,2/3,0];
Grpahene_TB_list.orbL = orbL;
%%
Grpahene_TB_list.show('HOPPING','scale', 2.4560000896,'atomscale',1,'TwoD',true);
%% Bulk band

t =1 ;
Grpahene_TB_n = Grpahene_TB.Subsall();
EIGENCAR = Grpahene_TB_n.EIGENCAR_gen();
bandplot(EIGENCAR);
%% 3d Band

Grpahene_TB_n.Rm = Rm;

%%
[EIGENCAR_3D,klist1,klist2] = Grpahene_TB_n.EIGENCAR_gen_3D([101,101],[-0.75,-0.25,0;1.5,-0.75,0;0 1.25 0]);
%%
[fig,ax] = Grpahene_TB_n.BZplot(Grpahene_TB_n.Rm,'mode','2D');
vasplib.bandplot_3d(EIGENCAR_3D,klist1,klist2,'ax',ax);