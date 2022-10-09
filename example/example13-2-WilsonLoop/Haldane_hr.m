clear;
sigma_0 = [1 0;0 1];sigma_x = [0 1;1 0];sigma_y = [0 -1i;1i 0];sigma_z = [1 0;0 -1];
Graphene = HR(2);
Graphene = Graphene<'POSCAR';
Graphene = Graphene<'KPOINTS';
search_range = [1 1 0];
maxR = 2.5;
Accuracy = 2;
Graphene = Graphene.nn(search_range, Accuracy ,maxR);
level_cut = 1;
Graphene = Graphene.H_TBSK_gen('level_cut',level_cut, ...
    'per_dir',[1 1 0]);
%%
syms m real;
Hal = Graphene.subs(sym('VppP_1'),sym('t','real'));
V_NNN = [-1 0 0; 0 -1 0;1 1 0];
Hal = Hal.set_hop(kron(sigma_z, -1i*m),1,1,V_NNN,'symadd');
Hal = Hal.autohermi();
%%
t = 1;
m = 0.1 * t;
Hal_n = Hal.Subsall();
EIGENCAR = Hal_n.EIGENCAR_gen();
bandplot(EIGENCAR ,[-3,3],'title',"Haldane-HR");
chern1 = chern_number(Hal_n,1);
chern2 = chern_number(Hal_n,2);
figure();
wilson_loop(Hal_n,"kx"); % This function may be wrongly set, raw edition from fyang 
%% vasplib inner function
chern0 = Hal_n.Chern()
chern1 = Hal_n.Chern('BAND_index',1)
chern2 = Hal_n.Chern('BAND_index',2)
%% check whether two results are same  
[BFCAR,~,klist_l] = Hal_n.WilsonLoop('knum_evol',101); % please check wanniertools answer
vasplib_plot.WilsonLoopPlot(BFCAR,klist_l)

[WannierCenterCAR,~,klist_l] = Hal_n.WannierCenter('knum_evol',101); % please check wanniertools answer
vasplib_plot.WccPlot(WannierCenterCAR,klist_l)