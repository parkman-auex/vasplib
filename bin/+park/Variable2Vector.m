function Vector = Variable2Vector(Variable)
%UNTITLED4 此处提供此函数的摘要
StrL = strsplit(string(Variable),'_');
StrL = StrL(2:4);
StrL = strrep(StrL,"m","-");
Vector = str2double(StrL);
end