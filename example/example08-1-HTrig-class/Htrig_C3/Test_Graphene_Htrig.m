%% Graphene
syms t a k_x k_y real;
bond = sqrt(3)/3 *a;
delta(1,:) = a*[1/2,-sqrt(3)/6];
delta(2,:) = a*[  0,    bond/a];
delta(3,:) = a*[-1/2,-sqrt(3)/6];
sigma_x = pauli_matric(1);
sigma_y = pauli_matric(2);
sigma_z = pauli_matric(3);
k = [k_x k_y].';
% Model
% 

Graphene = Htrig(2);
Graphene = Graphene +...
    Trig(-t*cos(delta(1,:)*k),sigma_x)+...
    Trig( t*sin(delta(1,:)*k),sigma_y)+...
    Trig(-t*cos(delta(2,:)*k),sigma_x)+...
    Trig( t*sin(delta(2,:)*k),sigma_y)+...
    Trig(-t*cos(delta(3,:)*k),sigma_x)+...
    Trig( t*sin(delta(3,:)*k),sigma_y);
Graphene_sym = Graphene.Htrig_sym
Graphene_sym = simplify(rewrite(Graphene_sym,'exp'))
% Band

a = 1;t =1;
Graphene_n  = Graphene.Subsall();
Graphene_n = Graphene_n<'POSCAR';
Graphene_n = Graphene_n <'KPOINTS';
EIGENCAR = Graphene_n.EIGENCAR_gen();
[klist_l,kpoints_l,kpoints_name] = Graphene_n.kpath_information();
bandplot(EIGENCAR,[-3,3],klist_l,kpoints_l,kpoints_name,'title','Graphene','Color','r');
%% Kane-Mele
%
% \begin{equation}
% \begin{aligned}
% &\begin{gathered}
% H=-t \sum_{\langle i, j\rangle \sigma} c_{\sigma \sigma}^{\dagger} c_{j \sigma}+i \lambda_{S O} \sum_{\langle\langle i, j\rangle\rangle \alpha, \beta} \nu_{i j} \sigma_{\alpha \beta}^{z} c_{i \alpha}^{\dagger} c_{j \beta} \\
% \quad+\lambda_{v} \sum_{i \sigma} \zeta_{i} c_{i \sigma}^{\dagger} c_{i \sigma} \cdot \\
% H=4 \lambda_{S O} \sin q_{x}\left(\cos q_{x}-\cos q_{y}\right) \Gamma^{15}-\lambda_{v} \Gamma^{2} \\
% \quad-t\left(2 \cos q_{x}+\cos q_{y}\right) \Gamma^{1}-t \sin q_{y} \Gamma^{12} \\
% \text { where } q_{x}=\frac{k_{x}}{2}, \text { and } q_{y}=\frac{\sqrt{3}}{2} k_{y} .
% \end{gathered}\\
% &\begin{array}{||c|c|c|c|c||}
% \hline \Gamma^{1} & \Gamma^{2} & \Gamma^{3} & \Gamma^{4} & \Gamma^{5} \\
% \tau^{x} \otimes I & \tau^{z} \otimes I & \tau^{y} \otimes \sigma^{x} & \tau^{y} \otimes \sigma^{y} & \tau^{y} \otimes \sigma^{z} \\
% \hline \multicolumn{5}{|l|}{\Gamma^{a b}=\left[\Gamma^{a}, \Gamma^{b}\right] /(2 i)} \\
% \hline
% \end{array}
% \end{aligned}
% \end{equation}
%% Model

syms t a lambda_SO lambda_v lambda_R Mz My Mx k_x k_y real;
syms q_x q_y;
sigma_z = pauli_matric(3);sigma_x = pauli_matric(1);
sigma_0 = pauli_matric(0);sigma_y = pauli_matric(2);
q_x = k_x/2;
q_y = sqrt(3)/2 *k_y;
%KaneMele = Htrig(4);
% try another type 
TYPE = 'sincos';% 'exp' 'sincos' 
KaneMele = Htrig(4,'Type',TYPE);
% gamma matrix maybe set wrongly!
KaneMele = KaneMele+...
    Trig(4*lambda_SO*sin(q_x)*(cos(q_x)-cos(q_y)),gamma_matric(1,5))+...
    Trig(-lambda_v                               ,gamma_matric(3  ))+...
    Trig(-t*(2*cos(q_x)+cos(q_y))                ,gamma_matric(1  ))+...
    Trig(-t*(sin(q_y))                           ,gamma_matric(1,3))+...
    Trig(4*lambda_SO*sin(q_x)*(cos(q_x)-cos(q_y)),gamma_matric(1,5))+...
    Trig(lambda_R*(1-cos(q_x)*cos(q_y))          ,gamma_matric(4  ))+...
    Trig(-sqrt(3)*lambda_R*sin(q_x)*sin(q_y)     ,gamma_matric(2  ))+...
    Trig(-lambda_R*cos(q_x)*sin(q_y)             ,gamma_matric(3,4))+...
    Trig(sqrt(3)*lambda_R*sin(q_x)*cos(q_y)      ,gamma_matric(3,2))+...
    Trig(Mz                                      ,sigma_z*sigma_0)+...
    Trig(Mx                                      ,sigma_x*sigma_0)+...
    Trig(My                                      ,sigma_y*sigma_0);

%% Band

a = 1;t =1;lambda_SO =0.02;lambda_v=0.0;lambda_R = 0.0;Mz = 0;Mx = 0;My = 0;
KaneMele_n  = KaneMele.Subsall();
KaneMele_n = KaneMele_n < 'POSCAR';
KaneMele_n = KaneMele_n <'KPOINTS';
EIGENCAR = KaneMele_n.EIGENCAR_gen();
[klist_l,kpoints_l,kpoints_name] = KaneMele_n.kpath_information();
bandplot(EIGENCAR,[-3,3],klist_l,kpoints_l,kpoints_name,'title','KaneMele','Color','r');
%% Descitize in C3
% KaneMele-zigzag
KaneMele.Rm = sym([1,0,0;-1/2,sqrt(3)/2,0;0 0 2/sqrt(3)])
KaneMele.orbL = [     0.333333333         0.666666667         0.50000000 ;...
    0.666666667         0.333333333         0.50000000;...
    0.333333333         0.666666667         0.50000000 ;...
    0.666666667         0.333333333         0.50000000...
    ];
KaneMele.Htrig_sym
KaneMele.rotation().Htrig_sym
KaneMele_dy = KaneMele.descritize([0,30,0],'Rotation','auto');
KaneMele_dy.Htrig_sym
%% 
a = 1;t =1;lambda_SO =0.06;lambda_v=0.4;lambda_R = -0.2;
KaneMele_dy = KaneMele_dy<'KPOINTS_slab';

[klist_l,kpoints_l,kpoints_name] = KaneMele_dy.kpath_information();
%KaneMele_dy.HcoeL = subs(KaneMele_dy.HcoeL);
% RCI_d = RCI_d.Subsall('para',[sym('theta_x'),sym('theta_y')]);
% [EIGENCAR_slab] = RCI_d.EIGENCAR_gen_slab('norb',-1,'paraname',[sym('theta_x'),sym('theta_y')],'para',[(0:0.1:1).',-(0:0.1:1).']);
KaneMele_dy_n = KaneMele_dy.Subsall();
[EIGENCAR_slab,WAVECAR_slab,WEIGHTCAR_slab] = KaneMele_dy_n.EIGENCAR_gen_slab();
%%
[fig,ax] = vasplib_tool.creat_figure();
vasplib.pbandplot(abs(WEIGHTCAR_slab),EIGENCAR_slab, ...
    'Ecut',[-2,2], ...
    'title','KaneMele-slab', ...
    'KPOINTS','KPOINTS_slab','fig',fig,'ax',ax);
%%
[fig,ax] = vasplib_tool.creat_figure();
vasplib.pbandplot(WEIGHTCAR_slab,EIGENCAR_slab, ...
    'Ecut',[-2,2], ...
    'title','KaneMele-slab', ...
    'KPOINTS','KPOINTS_slab', ...
    'cmap',ColorMap.Matplotlib('coolwarm'), ...
    'fig',fig,'ax',ax);
colorbar(ax);
%%
WaveFunc = WAVECAR_slab(:,59:60,60);
orb_list= KaneMele_dy.orbL;
%[fig,ax] = vasplib_tool.creat_figure();
vasplib.waveplot(orb_list,WaveFunc);

% Kane-Mele Armchair
% Dont know  how to solve them in only 4 orbital

KaneMele.Rm = sym([sqrt(3),1,0;-1/2,sqrt(3)/2,0;0 0 1/2])
KaneMele.orbL = [     0.333333333         0.666666667         0.50000000 ;...
    0.666666667         0.333333333         0.50000000;...
    0.333333333         0.666666667         0.50000000 ;...
    0.666666667         0.333333333         0.50000000...
    ];
KaneMele.Htrig_sym
B = KaneMele.rotation().Htrig_sym
%KaneMele_dy = KaneMele.descritize([0,30,0],'Rotation','auto');
%KaneMele_dy.Htrig_sym

%% 参考文献
%% 
% * <https://cpb-us-w2.wpmucdn.com/u.osu.edu/dist/3/67057/files/2018/09/graphene_tight-binding_model-1ny95f1.pdf 
% GrapheneTB(H(k) sin cos)>
% * <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.88.245115 KaneMeleTB-zigzag>
% * <https://arxiv.org/abs/1408.4507 KaneMeleTB-Armchair>
% * <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.48.11851 DescritizeTB 
% technique in square lattice>
% * Harper equation 