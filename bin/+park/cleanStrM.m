function StrM = cleanStrM(StrM)

cleanlabel = logical(1:size(StrM,1));
for i = 1:size(StrM,1)
    if ~isequal(StrM(i,:),repmat("",size(StrM(i,:))))
        cleanlabel(i) = false;
    end
end
StrM(cleanlabel,:)= [];
%
cleanlabel = logical(1:size(StrM,2));
for j = 1:size(StrM,2)
    if ~isequal(StrM(:,j),repmat("",size(StrM(:,j))))
        cleanlabel(j) = false;
    end
end
StrM(:,cleanlabel)= [];
end