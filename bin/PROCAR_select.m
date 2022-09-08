%% PROCAR_select
% usage WEIGHTCAR=PROCAR_select(PROCAR_collection,ionslist,orbitalslist,ionsnum,orbitalsnum)

function WEIGHTCAR=PROCAR_select(PROCAR_collection,ionslist,orbitalslist,bandsnum,kpointsnum)
%% to be extended to MX MY MZ
    
    WEIGHTCAR=zeros(bandsnum,kpointsnum);
    for i =1:length(ionslist)
        for j =1: length(orbitalslist)
            WEIGHTCAR=WEIGHTCAR+PROCAR_collection(ionslist(i),orbitalslist(j)).WEIGHTCAR;
        end
    end
end