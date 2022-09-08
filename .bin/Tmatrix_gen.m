function Tmatrix = Tmatrix_gen(H00,H01,w,eta)
    Dimi = length(H00);
    wc = w + 1i*eta;
    W = (wc*eye(Dimi)-H00);
    if abs(det(H01)) <1e-10
        H01 = H01 + 1e-6*eye(Dimi);
    end
    Tmatrix = [H01\W,-H01\(H01');eye(Dimi) zeros(Dimi)];
end