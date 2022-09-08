function after = tensor_symm(before,oper)
sizes = size(before);
rank = length(sizes);
after = before .* det(oper);
for i = 1:rank
    after = contract(after,1,oper,2);
end
end