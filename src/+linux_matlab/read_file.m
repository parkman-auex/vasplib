function [outstring] = read_file(filename)
import linux_matlab.*
    fid = fopen(filename, 'r');
    line_number = 0;
    outstring = [""];
    %-fgets 和 fgetl ： 可从文件读取信息
    while feof(fid) == 0
        line = fgetl(fid);
        outstring = [outstring ;line];
    end
    outstring(1,:) = [];
    fclose(fid);
end