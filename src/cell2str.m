function str_mat = cell2str(cell_mat)
[m, n] = size(cell_mat);
str_mat = strings(m,n);
for i = 1:m
    for j = 1:n
        str_mat(i,j) = string(cell_mat{i,j});
    end
end
end