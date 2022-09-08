function [Asort,Usort] = sorteig(U,A)
if nargin <2
    mode = 'eigenval';
    NUM_WAN = length(U);
else
    mode = 'whole';
    NUM_WAN = length(A);
end
NBANDS = length(U);

if strcmp(mode,'whole')
    % %--step9--:按从大到小的特征值顺序排序重新组合对应的特征向量
    Asort=zeros(NUM_WAN ,NBANDS );
    SortTmp=diag(U);%抽取特征值   
    [Usort,IJ]=sort(SortTmp,1,'ComparisonMethod','real');
    for jj=1:NBANDS 
        Asort(:,jj)=A(:,IJ(jj));%取特征向量的列向量        
    end
    Usort = diag(Usort);
elseif strcmp(mode,'eigenval')
    % %--step9--:按从大到小的特征值顺序排序重新组合对应的特征向量
    SortTmp=diag(U);%抽取特征值        
    [Usort,~]=sort(SortTmp,1);
    Usort = diag(Usort);
    Asort = [];
end