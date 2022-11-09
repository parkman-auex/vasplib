function Line = TotalLines(filename)
if ~isunix
    [status,OUTPUT] = system(strcat("find /V """" /C ",filename));
else
    Line = 0;
    warning('!USE IT IN WINDOWS');
end
if status
    Line = 0;
else
    OUTPUTCELL =  strsplit(OUTPUT,':');
    Line = str2double(OUTPUTCELL{2});
end
end