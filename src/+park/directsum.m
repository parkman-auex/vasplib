function C = directsum(A,B)
SizeA = size(A);
SizeB = size(B);
ZerosMat = zeros(SizeA(1),SizeB(2),class(A));
C = [A,ZerosMat ; ZerosMat.', B];
end