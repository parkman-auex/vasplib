classdef Term
    properties
        symbolic_polynomial;
        pauli_mat;
    end
    
    methods
        function term = Term(symbolic_polynomial,pauli_mat)
            if nargin < 1
                term.symbolic_polynomial = sym(0);
                term.pauli_mat = zeros(4);
            end
            %UNTITLED5 构造此类的实例
            %   此处显示详细说明
            term.symbolic_polynomial = symbolic_polynomial;
            term.pauli_mat = pauli_mat;
        end

        function term = disp(term)
            disp(sym(term));
        end
        function symmat = sym(term)
            termtmp = sym(zeros(length(term(1).pauli_mat.mat)));
            for i = 1:length(term)
                termtmp = termtmp + term(i).symbolic_polynomial * double(term(i).pauli_mat);
            end
            symmat = termtmp;
        end
        function C = plus(A,B)
            if isa(A,'Term') && isa(B,'Term')
                C = A;
                C(length(A)+1) = B;
            elseif isa(B,'Term') && ~isa(A,'Term')
                disp('gg');
                if isa(A,'HK')
                    C = A.plus(B);
                else
                    C = A;
                end
            elseif ~isa(B,'Term') && isa(A,'Term')
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
        function result = commute(A,B,isanti)
            if nargin < 3
                isanti = 0;
            end
            if isa(A,'pauli_matric') 
                A = A.mat;
            end
            if isa(B,'pauli_matric')
                B = B.mat;
            end
            M1 = sym(A);
            M2 = sym(B);
            if isanti == 1
                msymtmp = (M2*conj(M1))./((M1*M2));
            else
                msymtmp = (M2*M1)./(M1*M2);
            end
            mattmp = msymtmp/msymtmp(1,1);
            tf = diag(mattmp);
            label = all(tf == ones(length(tf),1));
            switch label
                case 1
                    eigentmp = diag(msymtmp);
                    result = eigentmp(1);
                case 0
                    disp('dependent')
                    if isanti == 1
                        result{1} = (M2*conj(M1))-(M1*M2);
                        result{2}  = (M2*conj(M1))+(M1*M2);
                    else
                        result{1} = (M2*M1)-(M1*M2);
                        result{2} = (M2*M1)+(M1*M2);
                    end
                    
            end
        end
        function term_list = subsk(term_list,kpoints_r)
            syms k_x k_y k_z;
            k_x = k_x - kpoints_r(1);
            k_y = k_y - kpoints_r(2);
            k_z = k_z - kpoints_r(3);
            for i = 1:length(term_list)
                term_list(i).symbolic_polynomial = subs(term_list(i).symbolic_polynomial);
            end
        end
    end
end

