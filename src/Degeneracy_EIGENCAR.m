function [DEGENCAR, NODEINSECT]= Degeneracy_EIGENCAR(EIGENCAR,highK,dE)
if nargin <3
    dE =0.1;
end
    [NBANDS,~] = size(EIGENCAR);
    NhighK = length(highK);
    DEGENCAR = zeros(NBANDS,NhighK);
    NODEINSECT = zeros(NBANDS-1,NhighK-1);
    % NhighK
    for i = highK
        DEGENCAR(1,i) = 1;
        count = 1;
        for j = 2:NBANDS
            if (EIGENCAR(j,i)-EIGENCAR(j-1,i)) < dE
                DEGENCAR(count,i) = DEGENCAR(count,i)+1;
            else
                count = count+1;
                DEGENCAR(count,i) = 1;
            end 
        end
    end
    % 
    for i = 1:NBANDS-1
        for j = 1:NhighK-1
            NODEINSECT(i,j) = sum(abs(EIGENCAR(i+1,highK(j)+1:highK(j+1)-1)-...
                EIGENCAR(i,highK(j):highK(j+1)-1)) < dE);
        end
    end
end