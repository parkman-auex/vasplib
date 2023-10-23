function Data = fileter1D(Data,windowsize)
    b = (1/windowsize)*ones(1,windowsize);
    a = 1;
    Data= filter(b,a,Data);
end