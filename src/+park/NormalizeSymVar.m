function SymVarList = NormalizeSymVar(SymVarList)
for i = 1:numel(SymVarList)
    A_SymVar = real(SymVarList(i));
    B_SymVar = imag(SymVarList(i));
    syms zzzz real ;
    A_SymVar_children = children(A_SymVar+zzzz);
    B_SymVar_children = children(B_SymVar+zzzz);
    A_SymVar_final = sym(0);B_SymVar_final = sym(0);
    for j = 1:length(A_SymVar_children)-1
        A_SymVar_final = A_SymVar_final + iNormalizeSymVar(A_SymVar_children{j});
    end
    for j = 1:length(B_SymVar_children)-1
        B_SymVar_final = B_SymVar_final + iNormalizeSymVar(B_SymVar_children{j});
    end
    SymVarList(i) = A_SymVar_final + B_SymVar_final*1i;
end
end

function SymVar = iNormalizeSymVar(SymVar,Accuracy)
if nargin <2
    Accuracy = 1e-8;
end

if isequal(simplify(SymVar),sym(0))
    SymVar = sym(0);
else
    [Coeffs,SvmVarBase] = coeffs(SymVar);   
    try
        CoeffsNum = double(Coeffs);
    catch 
        return;
    end
    % 
    SpecialBase = [1,sqrt(2),sqrt(3),sqrt(5)];
    SqecialPre = [1,2,3,4,5,6,7,8,9,...
        1/2,1/3,1/4,1/5,2/3,4/3,8/3,...
        ];
    for iSpecialBase = SpecialBase
        for jSqecialPre = SqecialPre
            CompareNum = jSqecialPre*iSpecialBase;
            if abs(CoeffsNum - CompareNum) < Accuracy
                SymVar = sym(CompareNum) * SvmVarBase;
                return;
            end
            if abs(CoeffsNum + CompareNum) < Accuracy
                SymVar = -sym(CompareNum) * SvmVarBase;
                return;
            end
        end
    end
end
end