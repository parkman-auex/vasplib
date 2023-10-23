function charL = num2strwithzero(numL)
charL = num2str(numL);
charL(charL == ' ') = '0';
end