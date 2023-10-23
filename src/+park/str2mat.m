function StrM = str2mat(StrL,delimiter,dim)
arguments
    StrL string;
    delimiter  = ' ';
    dim =2;
end
maxlength = 1;
nStrL = size(StrL,1);
for i = 1:nStrL
    %
    tmpStrL{i} = split(strtrim(StrL(i,:)),delimiter,dim);
    checklength{i} = length(tmpStrL{i});
    if checklength{i} > maxlength
        maxlength = checklength{i};
    else

    end
end
StrM = repmat("",[nStrL,maxlength]);
for i = 1:nStrL
    StrM(i,1:checklength{i}) = tmpStrL{i};
end
StrM = park.cleanStrM(StrM);
end