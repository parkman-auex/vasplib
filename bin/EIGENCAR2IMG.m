function IMG = EIGENCAR2IMG(EIGENCAR,dE,Erange,kselect)

% nargin  
if nargin < 2
    dE  = 0.1;
end
if nargin <3
    Erange = [min(min(EIGENCAR)),max(max(EIGENCAR))];
end
if nargin <4
    kselect = 1:size(EIGENCAR,2);
end
NBANDS = size(EIGENCAR,1);
DATA = EIGENCAR(:,kselect);

X_nodes = size(DATA,2);
Y_nodes = round(abs(Erange(2)-Erange(1))/dE);
Emin = min(Erange);
IMG_sparse = sparse(Y_nodes,X_nodes);
for i =1:X_nodes
    for j = 1:NBANDS
        NE = ceil((DATA(j,i)-Emin)/dE);
        if NE>0 && NE <= Y_nodes
            IMG_sparse(Y_nodes-NE+1,i) = IMG_sparse(Y_nodes-NE+1,i) +0.1;
        end
    end
end
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
IMG = full(IMG_sparse);
end

