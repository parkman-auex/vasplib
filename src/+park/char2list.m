function StrL = char2list(str)
arguments
    str char;
end
StrL = split(str,'\n');
checkL  = logical(1:numel(StrL));
for i = numel(StrL)
    if strcmp(StrL(i),'')
        checkL(i) = false;
    end
end
StrL = StrL(checkL);
end