clear
sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

Graphene = HR(2);
Graphene = Graphene<'POSCAR';
Graphene = Graphene<'KPOINTS';
search_range = [1 1 0];
maxR = 2.5;
Accuracy = 2;
Graphene = Graphene.nn_sk_smart(search_range, Accuracy ,maxR);
[Rnn,~,~,~] = Graphene.nn_information();
level_cut = 1;
Graphene = Graphene.H_TB_gen_SK('level_cut',1,'per_dir',[1 1 0]);
% list mode
Graphene_list = rewrite(Graphene);
Graphene = rewrite(Graphene_list,'rewind',true);
Graphene_list.symvar_list;
list(Graphene_list);
printout(Graphene);
%Graphene_list.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
Graphene_hk = Graphene_list.HR2HK();
sym(Graphene_hk);
%%
Graphene_htrig = Graphene_list.HR2Htrig();
latex(Graphene_htrig);
%%
K1 = [1/3,1/3,0];
Graphene_hk_K1 = Graphene_htrig.Htrig2HK(K1,'sym',true,'Order',1);
%%
K2 = [-1/3,-1/3,0];
Graphene_hk_K2 = Graphene_htrig.Htrig2HK(K2,'sym',true,'Order',1);
%%
G = [0,0,0];
Graphene_hk_G= Graphene_htrig.Htrig2HK(G,'sym',true,'Order',2);
%%
% C= simplify(sym(Graphene_htrig))
% simplify(C(2,1))
% 
% syms VppP_1 v_k k_x k_y k_plus real
% A= simplify(sym(Graphene_hk_K1))
% A = simplify(subs(A,[VppP_1] ,[2*v_k/(sqrt(3))] ))
% S = sym([exp(1i*pi/3) 0;0 exp(-1i*pi/3)])
% simplify(S*A*S' )
% 
% B = simplify(sym(Graphene_hk_K2))
% B = simplify(subs(B,[VppP_1] ,[2*v_k/(sqrt(3))] ))
% simplify(S'*B*S )
% 
% syms k phi_k real;
% A_prime =  simplify( rewrite(subs(A,[k_x k_y],[k*cos(phi_k) k*sin(phi_k)]),'exp'));
% A_prime = simplify(subs(A_prime,phi_k,phi_k+pi/3))
% 
% B_prime =  simplify( rewrite(subs(B,[k_x k_y],[k*cos(phi_k) k*sin(phi_k)]),'exp'));
% B_prime = simplify(subs(B_prime,phi_k,phi_k+pi/3))
% 
% D = subs(eig(S*A*S'),v_k,1);
% D(1) = -sqrt(D(1)^2);
% D(2) = sqrt(D(2)^2);
% fsurf(D)
% 
% syms VppP_1 v_k k_x k_y k_plus real
% H= simplify(sym(Graphene_hk_K3))
% H = simplify(subs(H,[VppP_1] ,[2*v_k/(sqrt(3))] ))
% S = sym([exp(1i*pi/6) 0;0 exp(-1i*pi/6)])
% simplify(S'*H*S)
