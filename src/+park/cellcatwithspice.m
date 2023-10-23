function OutputStr = cellcatwithspice(StrCell)
OutputStr = '';
if ischar(StrCell)
    OutputStr = strcat(OutputStr,StrCell,"\n");
    return;
end
for i = 1:length(StrCell)
    for j = 1:length(StrCell{i})-1
        OutputStr = strcat(OutputStr,StrCell{i}(j)," ");
    end
    OutputStr = strcat(OutputStr,StrCell{i}(j+1),"\n");
end
end