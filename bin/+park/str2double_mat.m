function [double_mat,str_mat] = str2double_mat(str)
    str = strtrim(str);
    str_list = strsplit(str,'\n');
    for i = 1:length(str_list)
        double_mat(i,:)=  str2double(strsplit(strtrim(str_list{i})));
    end
    % disp(str_mat);
    str_mat = string(double_mat);
end