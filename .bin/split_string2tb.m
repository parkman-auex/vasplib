%% toolfunction for split H_xyz

% syms k_1p k_12p k_2p k_3p k_1m k_12m k_2m k_3m a b c real

function [seq,vector,Hcoe,strout]=split_string2tb(strobj)
    vector =[0 0 0];
    containlist=["exp(k_1p)","exp(k_2p)","exp(k_3p)","exp(k_1m)","exp(k_2m)","exp(k_3m)","exp(k_12p)","exp(k_12m)","exp(k_1p2p)","exp(k_1p2m)"];

    if contains(strobj,containlist(1))
        vector=vector+[1,0,0];
        strobj=strrep(strobj,containlist(1),"1");
    end
    if contains(strobj,containlist(2))
        vector=vector+[0,1,0];
        strobj=strrep(strobj,containlist(2),"1");
    end
    if contains(strobj,containlist(3))
        vector=vector+[0,0,1];
        strobj=strrep(strobj,containlist(3),"1");
    end
    if contains(strobj,containlist(4))
        vector=vector+[-1,0,0];
        strobj=strrep(strobj,containlist(4),"1");
    end
    if contains(strobj,containlist(5))
        vector=vector+[0,-1,0];
        strobj=strrep(strobj,containlist(5),"1");
    end
    if contains(strobj,containlist(6))
        vector=vector+[0,0,-1];
        strobj=strrep(strobj,containlist(6),"1");
    end
    if contains(strobj,containlist(7))
        vector=vector+[1,0,0]+[0,-1,0];
        strobj=strrep(strobj,containlist(7),"1");
    end
    if contains(strobj,containlist(8))
        vector=vector+[-1,0,0]+[0,1,0];
        strobj=strrep(strobj,containlist(8),"1");
    end
    if contains(strobj,containlist(9))
        vector=vector+[1,0,0]+[0,1,0];
        strobj=strrep(strobj,containlist(9),"1");
    end
    if contains(strobj,containlist(10))
        vector=vector+[-1,0,0]+[0,-1,0];
        strobj=strrep(strobj,containlist(10),"1");
    end
    Hcoe=str2sym(strobj);
    v1=vector(1)+2;v2=vector(2)+2;v3=vector(3)+2;
    seq=(v1-1)*9+(v2-1)*3+v3;
    strout=strobj;
end

