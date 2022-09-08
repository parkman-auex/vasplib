function [outputArg1] = checkcell(input1,cellobj)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
outputArg1 = ~logical(1:length(cellobj));
for i = 1:length(cellobj)
    if  isequal(input1,cellobj{i})
        outputArg1 = true;
    end
end
end