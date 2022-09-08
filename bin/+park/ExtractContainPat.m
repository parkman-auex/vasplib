function [StrL1,StrL2] = ExtractContainPat(StrL,pat,options)
arguments
    StrL string{mustBeVector};
    pat = "=";
    options.IgnoreCase = false;
end

SelectL = logical(1:length(StrL));
for i =1:numel(StrL)
    if contains(StrL(i),pat,"IgnoreCase",options.IgnoreCase)
        SelectL(i) = false;
    end
end
StrL2 = StrL(SelectL);
StrL1 = StrL(~SelectL);
end