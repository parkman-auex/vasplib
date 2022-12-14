%% basis
% use 10.1103/PhysRevMaterials.2.054204 basis
w = exp(2*pi*1i/3);
% syms w;
psi{1} = (1i/sqrt(18))*[...
    1,1,-1,...
    w,w,-w,...
    w',w',-w',...
    1,1,-1,...
    w,w,-w,...
    w',w',-w'].';
psi{2} = (1i/sqrt(18))*[...
    1,1,-1,...
    w',w',-w',....
    w,w,-w,...
    1,1,-1,...
    w',w',-w',...
    w,w,-w].';
psi{3} = (1i/sqrt(18))*[...
    1,-1,-1,...
    -w,w,w,...
    w',-w',-w',...
    - 1,1,1,...
    w,-w,-w,...
    -w',w',w'].';
psi{4} = (1i/sqrt(18))*[...
    1,-1,-1,...
    - w',w',w',...
    w,-w,-w,...
    - 1,1,1,...
    w',-w',-w',...
    -w,w,w].';

%% test 
for i =1:4
  for j =1:4  
    fprintf('<psi%d|psi%d> = %4.2f\n',i,j,braket(psi{i},psi{j}));
  end
end
%%
phi{1} = (psi{1}+psi{3})/sqrt(2);
phi{2} = (psi{2}+psi{4})/sqrt(2);
phi{3} = (psi{1}-psi{3})/sqrt(2);
phi{4} = (psi{2}-psi{4})/sqrt(2);
%% +1 for A ;-1 for B
% use
% https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.123.256402/Supp_v0725.pdf
% sublattice
Subblattice_list  = [...
    1,-1,1,...
    -1,1,-1,...
    1,-1,1,...
    -1,1,-1,...s
    1,-1,1,...
    -1,1,-1 ...
    ];
Oper_S_psi = diag(Subblattice_list);

%% Oper_S for H
% Oper_S = <psi |Oper_S_psi|psi >
Oper_S = zeros(4);
for i = 1:4
    for j =1:4
        Oper_S(i,j) = braket(psi{i},Oper_S_psi,psi{j});
        fprintf('<psi%d|S|psi%d> = %4.2f\n',i,j,Oper_S(i,j));
    end
end
disp(Oper_S);
fprintf('*******************************\n The mat rep of S is\n');
disp(int8(real(Oper_S)));
%% Oper_S for H;
% Oper_S = <psi |Oper_S_psi|psi >
Oper_S2 = zeros(4);
for i = 1:4
    for j =1:4
        Oper_S2(i,j) = braket(phi{i},Oper_S_psi,phi{j});
        fprintf('<phi%d|S|phi%d> = %4.2f\n',i,j,Oper_S2(i,j));
    end
end
disp(Oper_S2);
fprintf('*******************************\n The mat rep of S is\n');
disp(int8(real(Oper_S2)));
%%
function value = braket(varargin)
    if nargin == 2
        psi1 = varargin{1};
        psi2 = varargin{2};
        value = psi1'*psi2;
    elseif nargin == 3
        psi1 = varargin{1};
        Oper = varargin{2};
        psi2 = varargin{3};
        value = psi1'*Oper*psi2;
    end
    
end