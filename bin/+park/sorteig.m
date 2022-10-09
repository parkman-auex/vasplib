function [Asort,Usort] = sorteig(U,A)
if nargin <2
    mode = 'eigenval';
    NUM_WAN = length(U);
else
    mode = 'whole';
    NUM_WAN = length(A);
end
NBANDS = length(U);
if ~isvector(U)
    SortTmp=diag(U);%抽取特征值
    vec = false;
else
    SortTmp = U;
    vec = true;
end
if strcmp(mode,'whole')
    % 按从大到小的特征值顺序排序重新组合对应的特征向量
    Asort=zeros(NUM_WAN ,NBANDS );
    [Usort,IJ]=sort(SortTmp,1,'ComparisonMethod','real');
    for jj=1:NBANDS
        Asort(:,jj)=A(:,IJ(jj));%取特征向量的列向量
    end
    if ~vec
        Usort = diag(Usort);
    end
elseif strcmp(mode,'eigenval')
    % 按从大到小的特征值顺序排序重新组合对应的特征向量
    SortTmp=diag(U);%抽取特征值
    [Usort,~]=sort(SortTmp,1,'ComparisonMethod','real');
    if ~vec
        Usort = diag(Usort);
    end
    Asort = [];
end
end