function [t_list,topo_of_Ham, topo] = UnsupClass_CM(...
    t_listA, topo_of_HamA, topoA,...
    t_listB, topo_of_HamB, topoB,...
    dim)
% clustering method of the results of function UnsupClass  
arguments
    t_listA double
    topo_of_HamA double
    topoA cell
    t_listB double
    topo_of_HamB double
    topoB cell
    dim int16
end
%% output1
t_list = [t_listA; t_listB];
%% output3
topo = topoA;
%% clustering method
LA = length(topoA);
LB = length(topoB);
for bi = 1:LB
    is_included = false;
    for ai = 1:LA        
        if isTopoEquv(topoA{ai},topoB{bi},dim)
            is_included = true;
            break
        end
    end
    
    if is_included
        topo_of_HamB = subs(topo_of_HamB,bi,-ai);
    else       
        topo(LA+1) = topoB(bi);
        topo_of_HamB = subs(topo_of_HamB,bi,-(LA+1));
    end
end
%% output2
topo_of_Ham = [topo_of_HamA, -topo_of_HamB];
end