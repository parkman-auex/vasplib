%% Retun 1 or 0 to show if there a string contains any string in string list

function flag = strcontain(str,containlist)

N=length(containlist);
flagtemp=0;
for i =1:N
    temp_add =  contains(str,string(containlist(i)));
    if strcmp(containlist(i),"")
        temp_add  = 1;
    end
    flagtemp = temp_add+flagtemp;
end

    flag =flagtemp;

end