classdef Trig
    properties
        symbolic_polynomial;
        pauli_mat;
    end
    
    methods
        function Trig = Trig(symbolic_polynomial,pauli_mat)
            if nargin < 1
                Trig.symbolic_polynomial = sym(0);
                Trig.pauli_mat = nan;
                Trig(1) = [];
            else
                Trig.symbolic_polynomial = symbolic_polynomial;
                Trig.pauli_mat = pauli_mat;
            end
        end
        function A = uminus(A)
            A.symbolic_polynomial = -A.symbolic_polynomial;
        end
        function C = plus(A,B)
            if isa(A,'Trig') && isa(B,'Trig')
                C = A;
                C(length(A)+1:length(A)+length(B)) = B;
            elseif isa(B,'Trig') && ~isa(A,'Trig')
                disp('gg');
                if isa(A,'HK')
                    C = A.plus(B);
                else
                    C = A;
                end
            elseif ~isa(B,'Trig') && isa(A,'Trig')
                if isa(B,'HK')
                    if length(B) == 2
                        C = B.plus(A);
                    else
                        C = B;
                    end
                else
                    C = B;
                end
            else
                C = 0;
            end
        end
        
    end
end

