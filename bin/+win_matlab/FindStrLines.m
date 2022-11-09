function Line = FindStrLines(filename,Str)
Str = strcat('"',Str,'"');
if ~isunix
    [status,OUTPUT] = system(strcat("find /N ",Str," ",filename));
else
    Line = 0;
    warning('!USE IT IN WINDOWS');
end
if status
    Line = 0;
else
    OUTPUTCELL =  strsplit(OUTPUT,{'[',']'});
    Line = str2double(OUTPUTCELL{2});
end
