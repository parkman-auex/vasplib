function [outstring,outnumber] = grep(filename, pattern,mode)
    %-模拟unix的grep指令
    %-filename：给出完整路径
    %-pattern:匹配表达式
    if nargin<3
        mode = 'aloud';
    end
    fid = fopen(filename, 'r');
    line_number = 0;
    outstring = [""];
    outnumber = [];
    %-fgets 和 fgetl ： 可从文件读取信息
    while feof(fid) == 0
        line = fgetl(fid);
        line_number = line_number + 1;
        matched = findstr(line, pattern);
        if ~isempty (matched)
            %-输出格式： 行号，对应行内容
            if strcmp(mode,'aloud')
                fprintf('%d: %s \n', line_number,line);
            end
            outstring = [outstring ;line];
            outnumber = [outnumber;line_number];
        end

    end
    outstring(1,:) = [];
    fclose(fid);
end