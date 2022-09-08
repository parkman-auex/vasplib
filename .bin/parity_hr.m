%%
% usage: Parity_list = Parity_hr(H_xyz,P_matrice,klist_s,mode)
%        Parity_list = Parity_hr(H_xyz,P_matrice,klist_s)
%        Parity_list = Parity_hr(H_xyz,P_matrice)
function Parity_list = Parity_hr(H_xyz,P_matrice,klist_s,mode,Accuracy)
%checkfirst
if nargin <5
Accuracy =4;
end
D_P = length(P_matrice);
if P_matrice * P_matrice ~= ones(D_P)
    error('PP != 1, run P matrice.');
end

if nargin<4
    mode = 't';
end

if nargin<3
    klist_s = [ 0 0 0;...
                0 0.5 0;...
                0.5 0 0;...
                0.5 0.5 0;
                0 0 0.5;...
                0 0.5 0.5;...
                0.5 0 0.5;...
                0.5 0.5 0.5];
end
[k_n,~] =size(klist_s);

Parity_list  =zeros(D_P,k_n);
if strcmp(mode,'t')
    % check value of orbital_init
    [norb ,~] =size(H_xyz(1).Hnum);
    if norb ~= D_P
        error("\n\n Oribital_of_P_matrice is wrong,please give a right P matrice!");
    end
    for i =1:k_n
        [~,WAVECAR] = EIGENCAR_gen(H_xyz,'m',0,-1,klist_s(i,:));
        [~,WAVECAR2] = EIGENCAR_gen(H_xyz,'m',0,-1,-klist_s(i,:));
        for j =1:D_P
            Parity_list(j,i) = OneOrZero(P_matrice,WAVECAR(:,j),WAVECAR2(:,j),Accuracy);
        end
    end
end

end
function [label] =  OneOrZero(P_matrice,vector1,vector2,Accuracy)
    
    vector3 = roundn(P_matrice*vector1,-Accuracy);
    vector2 = roundn(vector2,-Accuracy);
    if vector3 == vector2
        label = 1;
    elseif vector3 == -vector2
        label = -1;
    else
        label = 0;
        warning('something wrong?');
    end
end