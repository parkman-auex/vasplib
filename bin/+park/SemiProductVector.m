function KronVector = SemiProductVector(VectorLeft,VectorRight)
KronVector = [kron(VectorLeft,ones(size(VectorRight,1),1)),...
    repmat(VectorRight,[size(VectorLeft,1) 1])];
end