% Graphene
syms h theta phi real;
%
% phi   = 0;
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];
%%
syms a t delta k_x k_y real;
%
%
Rm = [1 0;-1/2 sqrt(3)/2]*a;
Graphene = HK(2,1);

t1 = [1/3,2/3]*a;
t2 = [2/3,1/3]*a;

e1 = t1-t2+[0  0]*Rm;
e2 = t1-t2+[1  0]*Rm;
e3 = t1-t2+[0 -1]*Rm;

% H12 = -t*(...
%     cos(e1*[k_x;k_y]) +cos(e2*[k_x;k_y])+cos(e3*[k_x;k_y])...
%     - ...
%     1i*sin(e1*[k_x;k_y]) +sin(e2*[k_x;k_y])+cos(e3*[k_x;k_y]));
% 
% H21 = conj(H12);
H12trig_r = Trig(-t*(cos(e1*[k_x;k_y]) +cos(e2*[k_x;k_y])+cos(e3*[k_x;k_y])),sigma_x);
H12trig_i = Trig(-t*(sin(e1*[k_x;k_y]) +sin(e2*[k_x;k_y])+sin(e3*[k_x;k_y])),sigma_y);
a = 1;
%
Graphene = Graphene + H12trig_r + H12trig_i+Term(delta,sigma_z) ;
Graphene.Rm = double(subs(Rm));
% [fig,ax] = BZplot(Graphene.Rm);

disp(Graphene);
%%
t = double(1);
delta = 0.5;
Graphene_n = Graphene.Subsall();
Graphene_n = Graphene_n < 'POSCAR';
Graphene_n = Graphene_n < 'KPOINTS';
%%
EIGENCAR = Graphene_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-3,3]);
%%
Graphene_hk = sym(Graphene_n);
Graphene_hk_n = subs(Graphene_hk);
[W,E] = eig(Graphene_hk_n);
%%
Bcsym = vasplib.Berry_curvature_kubo((Graphene_hk_n),k_x,k_y);
Bcfun = matlabFunction(Bcsym,'Var',[k_x,k_y]);
%%
count = 0;
size = 101;
R = zeros(size*size,3);
A = zeros(size*size,3);
for kx = linspace(-2*pi,2*pi,size)
    for ky = linspace(-2*pi,2*pi,size)
        count =count +1;
        R(count, :) = [ kx,ky,0];
        Bc_n = Bcfun(kx,ky);
        A(count,:) = [0,0,Bc_n(2)];
    end
end
[fig,ax] = BZplot(Graphene.Rm);
A(:,3) = A(:,3)*1e20;
vasplib_tool.quiverplot(R,A,'b','Graphene',fig,ax)
% klist = rand(4,2);
% Bc_n = Berry_curvature_2bandn(subs(Graphene_kp),k_x,k_y,klist)
% u1 = simplify(NomalizeEigenvector(W(:,2))); 

% Bc = Berry_curvature_2D(u1,k_x,k_y);
%%
% the defination way is really difficult
% use anotherway for two level system as a test


