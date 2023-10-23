classdef quaternion
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        U

    end
    properties(Dependent)
        expression
    end

    methods
        function Q = quaternion(Re,I,J,K)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            Q.U = [Re;I;J;K];
        end
    end
    % get
    methods
        function expression = get.expression(Q)
            syms I J K;
            expression = Q.U(1)+Q.U(2)*I+Q.U(3)*J+Q.U(4)*K;
        end
    end
    methods(Static)
        function Q_I = I()
            Q_I = quaternion(0,1,0,0);
        end
        function Q_I = J()
            Q_I = quaternion(0,0,1,0);
        end
        function Q_I = K()
            Q_I = quaternion(0,0,0,1);
        end
    end
    methods
        function Q = conj(Q)
            Q.U(2:4) = -Q.U(2:4);
        end
        function Value = innerproduct(Q1,Q2)
            Value = norm(Q1*conj(Q2));
        end
        function normvalue = norm(Q)
            normvalue = norm(Q.U);
        end
        function Q =  inverse(Q)
            Q = 1/norm(Q).^2 * conj(Q);
        end
        function Re = real(Q)
            Re = Q.U(1);
        end
        function Q = mtimes(Q1,Q2)
            % \begin{aligned}
            % \hat{\mathbf{q}} \hat{\mathbf{r}} & =\left(i q_x+j q_y+k q_z+q_w\right)\left(i r_x+j r_y+k r_z+r_w\right) \\
            % & =i\left(q_y r_z-q_z r_y+r_w q_x+q_w r_x\right) \\
            % & +j\left(q_z r_x-q_x r_z+r_w q_y+q_w r_y\right) \\
            % & +k\left(q_x r_y-q_y r_x+r_w q_z+q_w r_z\right) \\
            % & +q_w r_w-q_x r_x-q_y r_y-q_z r_z \\
            % & =\left(\hat{\mathbf{q}}_v \times \hat{\mathbf{r}}_v+r_w \hat{\mathbf{q}}_v+q_w \hat{\mathbf{r}}_v, q_w r_w-\hat{\mathbf{q}}_v \cdot \hat{\mathbf{r}}_v\right)
            % \end{aligned}
            if isa(Q1,'quaternion') && isa(Q2,'quaternion')
                Q = Q1;
                q_x = Q1.U(2);
                q_y = Q1.U(3);
                q_z = Q1.U(4);
                q_w = Q1.U(1);

                r_x = Q2.U(2);
                r_y = Q2.U(3);
                r_z = Q2.U(4);
                r_w = Q2.U(1);
                Q.U(1) = q_w * r_w-q_x * r_x-q_y * r_y-q_z * r_z;
                Q.U(2) = q_y * r_z-q_z * r_y+r_w * q_x+q_w * r_x;
                Q.U(3) = q_z * r_x-q_x * r_z+r_w * q_y+q_w * r_y;
                Q.U(4) = q_x * r_y-q_y * r_x+r_w * q_z+q_w * r_z;

            elseif isa(Q1,'quaternion') || ~isa(Q2,'quaternion')
                Q = Q1;
                Q.U = Q.u * Q2;

            elseif ~isa(Q1,'quaternion') && isa(Q2,'quaternion')
                Q = Q2;
                Q.U = Q.u * Q1;

            end
        end
        function Q = plus(Q1,Q2)
            Q= Q1;
            Q.U = Q.U+Q2.U;
        end
        function Q = uminus(Q)
            Q.U = -Q.U;
        end
        function Q = minus(Q1,Q2)
            Q = Q1 + (-Q2);
        end
        function expressionform  =disp(Q)
            expressionform  = Q.expression;
            disp(expressionform);
        end
    end
end