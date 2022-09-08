function StrM = char2mat(str)
arguments
    str char;
end
StrL = split(str,'\n');
StrM = park.str2mat(StrL);
end