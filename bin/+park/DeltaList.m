function OutPutM = DeltaList(LeftL,RightL)
    % Width = size(LeftL,2);
    Ncol = size(RightL,1);
    OutPutM = zeros(Ncol);
    [ia,ib] = ismember(LeftL,RightL,'rows');
    iL = find(ia);
    jL = ib(ia);
    for i = 1:length(iL)
        OutPutM(iL(i),jL(i)) = 1;
    end
end