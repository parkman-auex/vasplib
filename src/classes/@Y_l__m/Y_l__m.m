classdef Y_l__m < HollowKnight
    %Ylm  Spherical harmonic function.
    %   Spherical harmonics of order l (0<=l<=bandlimit) and degree m (-l<=m<=l):
    %
    %                                         ( 2l+1   (l-m)! )^(1/2)
    %       Y_{l}^{m} (theta, phi) = (-1)^m * ( ---- * ------ )       * exp(i*m*phi) * P_{l}^{m} (cos(theta))
    %                                         ( 4*pi   (l+m)! )
    %
    %   where P_{l}^{m} (x) is the associated Legendre polynomial of order l
    %   and degree m.
    %   Y = YLM(L,M) computes the  spherical harmonic of orbital angular momentum
    %   L and magnetic angular momentum  M
    %
    %   Y = YLM(L,M,TH,PHI) computes the surface spherical harmonic of
    %   degree L and order M, evaluated at each element of inclination TH and
    %   azimuth PHI. N and M must be scalar integers where M <= abs(N). TH and
    %   PHI must be arrays of equal size.
    %
    %
    %   Y = YLM(__,Name,Value) specifies additional options using one or
    %   more name-value pair arguments. Valid options are:
    %
    %     - 'symbolic' specifies whether to compute the complex spherical harmonics
    %       or their real part. Valid values are 'complex' (default) and
    %       'real'. Real spherical harmonics are of cosine type for M > 0 and
    %       of sine type for M < 0.
    %
    %     - 'normalization' specifies whether the result of the computation is to be
    %       normalized. The normalization coefficient is chosen so as to ensure
    %       that the spherical harmonics are orthonormal. Valid values are true
    %       (default), false.
    %
    %     - 'phase' specifies whether to include the Condon-Shortley phase
    %       term. This term is not strictly necessary but may be useful in
    %       quantum mechanics applications. Valid values are true (default),
    %       false.
    %
    %   See also LEGENDRE.
    properties
        l   double;
        m   double;
        % coe ;
    end
    properties(Dependent)
        expression ;
        parity ;
    end
    properties(Dependent,Hidden)
        hollow;
    end
    properties(Hidden)
        n   = 0         ;
        symbolic = true ;
        explicit = false;
    end
    methods
        function YlmObj = Y_l__m(l,m,coe,n,TH,PHI,options)
            arguments
                l   = 1;%... Azimuthal quantum number
                m   = [];%... Magnetic quantum number
                coe = 1;
                n   = 0;%... Radius(n?)
                TH  = sym('theta','real');%... azimutal angles
                PHI = sym('phi','real');%... polar angles
                options.common = true;
                options.symbolic = false;
                options.normalization {mustBeMember(options.normalization,{'unnorm','norm','sch'})}= 'unnorm';
            end

            % call father : HollowKnight
            YlmObj = YlmObj@HollowKnight();
            
            % single
            YlmObj.symbolic = options.symbolic ;
            if ~options.common
                if isa(l,'sym') && isa(m,'sym')
                    switch options.normalization
                        case 'unnorm'
                            YlmObj.expression = sym('');
                        case 'norm'
                            YlmObj.expression = sym('');
                        case 'sch'
                    end
                else
                    
                end
                if isa(TH,'sym')
                    
                end
                if isa(PHI,'sym')
                end
            end
            if isempty(m)
                m = (l:-1:-l).';
            end
            YlmObj(1,1).l = l(1,1);
            YlmObj(1,1).m = m(1,1);
            YlmObj(1,1).coe = coe(1,1);
            YlmObj(1,1).n = n(1,1);

            % multi
            [nS_row,nS_col] = size(l);
            [nSz_row,nSz_col] = size(m);
            [ncoe_row,ncoe_col] = size(coe);
            [nn,nn_col] = size(n);
            nYlmObj = max([nS_row ,nSz_row, ncoe_row,nn]);
            nYlmObj_col = max([nS_col ,nSz_col, ncoe_col nn_col] );
            if nYlmObj == 1 && nYlmObj_col == 1
                return;
            else
                YlmObj = repmat(YlmObj(1,1),[nYlmObj nYlmObj_col]);
                [l,m,coe,n] = Y_l__m.StanderdInput(l,m,coe,n,'nrow',nYlmObj,'ncol',nYlmObj_col);
                for i =1:numel(YlmObj)
                    YlmObj(i).l = l(i);
                    YlmObj(i).m = m(i);
                    YlmObj(i).n = n(i);
                    YlmObj(i).coe = coe(i);
                end  
            end
        end
    end
    
    methods %get
        function hollow = get.hollow(YlmObj)
            if isnan(YlmObj.coe)
                hollow = true;
            else
                hollow = false;
            end
        end
        function parity = get.parity(YlmObj)
            if length(YlmObj) == 1
                parity = (-1)^YlmObj.l;
            else
                parity = YlmObj(1).parity;
                for i =2:length(YlmObj)
                    if parity ~=YlmObj(i).parity
                        parity = 0 ;
                        return;
                    end
                end
            end
        end
        function expression = get.expression(YlmObj)
            if YlmObj.l>4
                syms theta phi real;
                if YlmObj.symbolic
                    Symble = vasplib.SymbolicVarible('N',YlmObj.m,YlmObj.l);
                    expression = Symble...
                        *((-1)^YlmObj.m*(sqrt(1/2/sym(pi)))*exp(1i*YlmObj.m*phi));
                    expression = simplify(YlmObj.coe*expression);
                else
                    Symble = vasplib.SymbolicVarible('N',YlmObj.m,YlmObj.l);
                    expression = Symble...
                        *vpa((-1)^YlmObj.m*sqrt(1/2/pi)*exp(1i*YlmObj.m*phi));
                end
            else
                expression = explicitformula(YlmObj,'vpa',false);
                expression = simplify(YlmObj.coe*expression);
            end
        end
    end
    methods % reload
        function C = plus(A,B)
            if isa(A,'Y_l__m') && isa(B,'Y_l__m')
                C = A;
                C = [C,B];
                C = contract(C);
            end
        end
        function A = umius(A)
            if isa(A,Y_l__m)
                for i =1:length(A)
                    A(i).coe = -A(i).coe;
                end
            end
        end
        function C = minus(A,B)
            C = A + -B;
        end
        function C = innertimes(A,B)
            if A.n ~= B.n
                C = 0;
                return;
            end
            l1 = A.l;l2 = B.l;
            m1 = A.m;m2 = B.m;
            % sqrt{\frac{(2 a+1)(2 b+1)}{4 \pi}}
            % selection rule m1 + m2 = - m3
            M = m1+m2;
            count = 0;
            if A.symbolic ||B.symbolic
                CoePublic = sqrt((2*l1+1)*(2*l2+1)/sym(pi))/2 ;
                for L = abs(l1-l2):abs(l1+l2)
                    C1 = sym((-1)^L*sqrt(2*L+1));
                    W1 = Y_l__m.Wigner3j_sym([l1 l2 L],[m1 m2 -M]);
                    W2 = Y_l__m.Wigner3j_sym([l1 l2 L],[0 0 0]);
                    if W1~=sym(0) && W2~=sym(0)
                        count = count +1;
                        lL(count) = L;
                        COE(count) = C1*W1*W2;
                    end
                    %m(count) = m;
                end
            else
                CoePublic = sqrt((2*l1+1)*(2*l2+1)/pi)/2 ;
                for L = abs(l1-l2):abs(l1+l2)
                    C1 = (-1)^L*sqrt(2*L+1);
                    W1 = Y_l__m.Wigner3j([l1 l2 L],[m1 m2 -M]);
                    W2 = Y_l__m.Wigner3j([l1 l2 L],[0 0 0]);
                    if W1~=0 && W2~=0
                        count = count +1;
                        lL(count) = L;
                        COE(count) = C1*W1*W2;
                    end
                    %m(count) = m;
                end
            end
            C  = Y_l__m(lL,M,CoePublic*COE,A.n);

        end
        function C = eq(A,B,options)
            arguments
                A Y_l__m;
                B Y_l__m;
                options.strict = false;
            end
            C = true;
            if length(A) == 1 && length(B) == 1
                if A.l ~= B.l
                    C = false;
                    return;
                end
                if A.m ~= B.m
                    C = false;
                    return;
                end
                if A.n ~= B.n
                    C = false;
                    return;
                end
                if options.strict
                    if A.coe ~= B.coe
                        C = false;
                        return;
                    end
                end
            elseif isvector(A) && isvector(B)
                A = contract(A);
                B = contract(B);
                if length(A) ~= length(B)
                    C = false;return;
                end
                if isequal([A(:).l],[B(:).l])
                    C = false;return;
                end
                if isequal([A(:).m],[B(:).m])
                    C = false;return;
                end
                if isequal([A(:).n],[B(:).n])
                    C = false;return;
                end
                % linearly independent rule
                CoeLtmp = [A(:).coe;B(:).coe];
                if rank(CoeLtmp) == 2
                    C = false;return;
                end
                if options.strict
                    if isequal(CoeLtmp(1,:),CoeLtmp(2,:))
                        C = false;
                        return;
                    end
                end
            else
                C = false;
            end
        end
        function C = mrdivide(A,B)
            A = contract(A);
            B = contract(B);
            if A == B
                C = A(1).coe/B(1).coe;
            else
                C = 0;
            end
        end
        function C = lt(A,B,options)
            arguments
                A Y_l__m;
                B Y_l__m;
                options.strict = false;
            end
            C = true;
            if length(A) == 1 && length(B) == 1
                if A.n > B.n
                    C = false;return;
                end
                if A.l > B.l
                    C = false;return;
                end
                if A.m > B.m
                    C = false;return;
                end
                if options.strict
                    if A == B
                        if A.coe >= B.coe
                            C = false;return;
                        end
                    end
                end
                %return;
            else
                % need to check
                C = length(A) < length(B);
            end
        end
        function sum()
            
        end
    end
    %% contract
    methods 
        function B = contractrow(A,options)
            arguments
                A Y_l__m;
                options.forgetcoe = false;
            end
            %
            B = A;%keep origin
            %cleanzero first
            A = cleanrow(A);
            % delete zero coe!!
            lmL = ([A(:).l;A(:).m;]).';
            coeL = [A(:).coe];
            [lmL_unique,~,~] = unique(lmL,'rows');
            if isempty(lmL_unique)
                %B = A;
                return;
            end
            if options.forgetcoe
                A = repmat(A(1),[1 length(lmL_unique)]);
                for i =1:length(lmL_unique)
                    A(i).l = lmL_unique(i,1);
                    A(i).m = lmL_unique(i,2);
                    A(i).coe = 1;
                end
            else
                B = repmat(A(1),[1 size(lmL_unique,1)]);
                for i = 1:size(lmL_unique,1)
                    lml_tmp = lmL_unique(i,:);
                    [~,index_ic] = ismember(lmL,lml_tmp,'rows');
                    B(1,i).l = lmL_unique(i,1);
                    B(1,i).m = lmL_unique(i,2);
                    B(1,i).coe = sum(coeL(logical(index_ic)));
                end
            end
            B = cleanrow(B);
        end
    end
    %% rotation
    methods
        function Ak = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
            arguments
                A
                abc
                RightorLeft
                immproper = false;
                conjugate = false;
                antisymmetry = false;
            end
            alpha = (abc(1));
            beta = (abc(2));
            gamma = (abc(3));
            A_L = Y_l__m(A.l);
            Ak = A.HollowMe;
            for i =  1:length(A_L)
                Ai = A_L(i);
                if A.l == 0
                    Ai.coe = 1;
                else
                    m1 = A.m;
                    m2 = Ai.m;
                    WignerD_single_element = (Y_l__m.d(A.l,m1,m2,beta));
                    Ai.coe = Ai.coe*exp(1i*RightorLeft*m1*alpha)*WignerD_single_element*exp(1i*RightorLeft*m2*gamma);
                    % for Y_l__m
                    Ai.coe = conj(Ai.coe);
                end
                if immproper
                    Ai.coe = (-1)^(A.l)*Ai.coe;
                end
                if ~conjugate
                    Ai.coe = Ai.coe * A.coe;
                else
                    Ai.coe = conj(Ai.coe * A.coe);
                end
                Ak = [Ak,Ai];
            end
        end
    end
    %% others
    methods % disp
        function disp(YlmObj,options)
            arguments
                YlmObj Y_l__m;
                options.vpa = true;
                options.explicit = true;
                options.cart = false
            end
            optionsCell = namedargs2cell(options);
            if length(YlmObj)== 1
                disp(YlmObj.explicitformula(optionsCell{:}));
            else
                disp(YlmObj.formula(optionsCell{:}));
            end
            
        end
        function expression = pretty(YlmObj,options)
            
        end
        function Tesseralexpansion = Tesseral(Y_lmObj,options)
            arguments
                Y_lmObj;
                options.sym = true;
            end
            % take care! Y_lm is actually Y_l__m with a vectorize
            % combination
            % we use Tesseral get real SH coe
            % Y_{\ell}^{m}= \begin{cases}\frac{1}{\sqrt{2}}\left(Y_{\ell|m|}-i Y_{\ell,-|m|}\right) & \text { if } m<0 \\ Y_{\ell 0} & \text { if } m=0 \\ \frac{(-1)^{m}}{\sqrt{2}}\left(Y_{\ell|m|}+i Y_{\ell,-|m|}\right) & \text { if } m>0\end{cases}
            for i = 1:size(Y_lmObj,1)       % l m coe
                Y_lmObjtmp = Y_lmObj(i,:);
                comparerowL = [];
                if options.sym
                    sumrowL = sym([]);
                else
                    sumrowL = [];
                end
                for j = 1:length(Y_lmObjtmp)
                    l = Y_lmObjtmp(j).l;
                    m = Y_lmObjtmp(j).m;
                    %
                    if m < 0
                        coe1 = 1/sqrt(2);
                        coe2= -1i/sqrt(2);
                    elseif m>0
                        coe1 = 1/sqrt(2)*(-1)^m;
                        coe2 = 1i/sqrt(2)*(-1)^m;
                    else
                        coe1 = 1/2;
                        coe2 = 1/2;
                    end
                    comparerowL = [comparerowL;l abs(m);l -abs(m)];
                    sumrowL = [sumrowL;coe1*Y_lmObjtmp(j).coe;coe2*Y_lmObjtmp(j).coe];
                end
                [Tesseralexpansion{i,1},Tesseralexpansion{i,2}] = Y_l__m.generalcontractrow(comparerowL,sumrowL);
            end
        end
        function str = string(YlmObj,options)
            arguments
                YlmObj Y_l__m;
                options.vpa = true;
                options.explicit = false;
                options.cart = true;
            end
            optionsCell = namedargs2cell(options);
            str = string(formula(YlmObj,optionsCell{:}));
        end
        function expression = formula(YlmObj,options)
            arguments
                YlmObj Y_l__m;
                options.vpa = true;
                options.explicit = false;
                options.cart = false;
            end
            expression =sym(0);
            optionsCell = namedargs2cell(options);
            if options.cart
                expression = repmat(expression,[size(YlmObj,1),1]);
                for j = 1:size(YlmObj,1)
                    %expressiontmp =sym(0);
                    for i = 1:size(YlmObj,2)
                        if ~YlmObj(j,i).hollow
                            expression(j) = expression(j)  + YlmObj(j,i).coe*Y_l__m.nlm2atomic(YlmObj(j,i).l,YlmObj(j,i).m,YlmObj(j,i).n,'outputformat','sym');
                        end
                    end
                    %expression = [expression;expressiontmp];
                end
                expression = simplify(expression);
                return;
            end
            if options.explicit
                for i = 1:length(YlmObj)
                    expression = expression +YlmObj(i).coe*explicitformula(YlmObj(i),optionsCell{:});
                end
            else
                if options.vpa
                    for i = 1:length(YlmObj)
                        Symble = vasplib.SymbolicVarible('Y',YlmObj(i).m,YlmObj(i).l);
                        expression = expression + vpa(YlmObj(i).coe)*Symble;
                    end
                else
                    for i = 1:length(YlmObj)
                        Symble = vasplib.SymbolicVarible('Y',YlmObj(i).m,YlmObj(i).l);
                        expression = expression + YlmObj(i).coe*Symble;
                    end
                end
            end
        end
        function expression = explicitformula(YlmObj,options)
            arguments
                YlmObj Y_l__m;
                options.vpa = true;
                options.explicit = true;
                options.cart = false
                options.seed =["theta","phi"];
            end
            optionsCell = namedargs2cell(options);
            if YlmObj.m < 0
                %YlmObj2 = YlmObj;
                YlmObj.m = -YlmObj.m;
                % P_{\ell }^{-m}(x)=(-1)^{m}{\frac {(\ell -m)!}{(\ell +m)!}}P_{\ell }^{m}(x).
                %expression = -1^YlmObj.m*factorial(YlmObj.l-YlmObj.m)/factorial(YlmObj.l+YlmObj.m)...
                %    *explicitformula(YlmObj,options{:});
                expression = conj(explicitformula(YlmObj,optionsCell{:}))*(-1)^YlmObj.m;
                return;
            end
            if ~options.cart
                theta = sym(options.seed(1),'real');
                phi = sym(options.seed(2),'real');
            else
                syms x y z r real;
            end
            % <https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#cite_note-Varshalovich1988-1>
            if options.cart
                coe_pi = sqrt(1/sym(pi));
                coe_m = (-1)^YlmObj.m;
                coe_r = 1/r^YlmObj.l;
                XIYfunc = (x+1i*y)^YlmObj.m;
                switch  YlmObj.l
                    case 0
                        coe_1 = 1/2;%1/2
                        coe_front_pi = 1;
                        Zfunc = 1;
                    case 1
                        coe_1 = 1/2;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(3);
                                Zfunc = z; 
                            case 1
                                coe_front_pi = sqrt(3/2);
                                Zfunc = 1;
                        end
                    case 2
                        coe_1 = 1/4;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(5);
                                Zfunc = 3*z^2-r^2;
                            case 1
                                coe_front_pi = sqrt(30);
                                Zfunc = z; 
                            case 2
                                coe_front_pi = sqrt(15/2);
                                Zfunc = 1;
                        end
                    case 3
                        coe_1 = 1/4;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(7);
                                Zfunc = (5*z^2-3*r^2)*z;
                            case 1
                                coe_front_pi = sqrt(21/4);
                                Zfunc = (5*z^2-r^2);
                            case 2
                                coe_front_pi = sqrt(105/2);
                                Zfunc = z; 
                            case 3
                                coe_front_pi = sqrt(35/4);
                                Zfunc = 1;
                        end
                    case 4
                        coe_1 = 3/16;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(1);
                                Zfunc = (35*z^4-30*z^2*r^2+3*r^2);
                            case 1
                                coe_front_pi = sqrt(40);
                                Zfunc = (7*z^2-3*r^2)*z;
                            case 2
                                coe_front_pi = sqrt(8);
                                Zfunc = (7*z^2-r^2);
                            case 3
                                coe_front_pi = sqrt(35*4);
                                Zfunc = z; 
                            case 4
                                coe_front_pi = sqrt(35/2);
                                Zfunc = 1;
                        end
                    case 5
                        coe_1 = 1/16;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(11);
                                Tfunc = SINT^YlmObj.m*(63*COST^5-70*COST^3+15*COST);
                            case 1
                                coe_front_pi = sqrt(165/2);
                                Tfunc = SINT^YlmObj.m*(21*COST^4-14*COST^2+1);
                            case 2
                                coe_front_pi = sqrt(1155*2);
                                Tfunc = SINT^YlmObj.m*(3*COST^3-COST);
                            case 3
                                coe_front_pi = sqrt(385/4);
                                Tfunc = SINT^YlmObj.m*(9*COST^2-1);
                            case 4
                                coe_front_pi = sqrt(385*9/2);
                                Tfunc = SINT^YlmObj.m*COST;
                            case 5
                                coe_front_pi = sqrt(77*9/4);
                                Tfunc = SINT^YlmObj.m;
                        end
                    otherwise
                        expression =Y_l__m.InnerY_l__m(YlmObj.l,YlmObj.m,'triangle',false);
                        return;
                        
                end
                expression = coe_m*coe_1*coe_front_pi*coe_pi*XIYfunc*Zfunc*coe_r;
            else
                EIPHI = exp(1i*YlmObj.m*phi);
                COST = cos(theta);SINT = sin(theta);
                coe_pi = sqrt(1/sym(pi));
                coe_m = (-1)^YlmObj.m;
                switch  YlmObj.l
                    case 0
                        coe_1 = 1/2;%1/2
                        coe_front_pi = 1;
                        Tfunc = 1;
                    case 1
                        coe_1 = 1/2;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(3);
                                Tfunc = COST;
                            case 1
                                coe_front_pi = sqrt(3/2);
                                Tfunc = SINT;
                        end
                    case 2
                        coe_1 = 1/4;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(5);
                                Tfunc = 3*COST^2-1;
                            case 1
                                coe_front_pi = sqrt(30);
                                Tfunc = SINT*COST;
                            case 2
                                coe_front_pi = sqrt(15/2);
                                Tfunc = SINT^2;
                        end
                    case 3
                        coe_1 = 1/4;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(7);
                                Tfunc = SINT^YlmObj.m*(5*COST^3-3*COST);
                            case 1
                                coe_front_pi = sqrt(21/4);
                                Tfunc = SINT^YlmObj.m*(5*COST^2-1);
                            case 2
                                coe_front_pi = sqrt(105/2);
                                Tfunc = SINT^YlmObj.m*COST;
                            case 3
                                coe_front_pi = sqrt(35/4);
                                Tfunc = SINT^YlmObj.m;
                        end
                    case 4
                        coe_1 = 3/16;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(1);
                                Tfunc = SINT^YlmObj.m*(35*COST^4-30*COST^2+3);
                            case 1
                                coe_front_pi = sqrt(40);
                                Tfunc = SINT^YlmObj.m*(7*COST^3-3*COST);
                            case 2
                                coe_front_pi = sqrt(8);
                                Tfunc = SINT^YlmObj.m*(7*COST^2-1);
                            case 3
                                coe_front_pi = sqrt(35*4);
                                Tfunc = SINT^YlmObj.m*COST;
                            case 4
                                coe_front_pi = sqrt(35/2);
                                Tfunc = SINT^YlmObj.m;
                        end
                    case 5
                        coe_1 = 1/16;%1/2
                        switch YlmObj.m
                            case 0
                                coe_front_pi = sqrt(11);
                                Tfunc = SINT^YlmObj.m*(63*COST^5-70*COST^3+15*COST);
                            case 1
                                coe_front_pi = sqrt(165/2);
                                Tfunc = SINT^YlmObj.m*(21*COST^4-14*COST^2+1);
                            case 2
                                coe_front_pi = sqrt(1155*2);
                                Tfunc = SINT^YlmObj.m*(3*COST^3-COST);
                            case 3
                                coe_front_pi = sqrt(385/4);
                                Tfunc = SINT^YlmObj.m*(9*COST^2-1);
                            case 4
                                coe_front_pi = sqrt(385*9/2);
                                Tfunc = SINT^YlmObj.m*COST;
                            case 5
                                coe_front_pi = sqrt(77*9/4);
                                Tfunc = SINT^YlmObj.m;
                        end
                    otherwise
                        expression =Y_l__m.InnerY_l__m(YlmObj.l,YlmObj.m);
                        return;
                end
                expression = coe_m*coe_1*coe_front_pi*coe_pi*Tfunc*EIPHI;
            end
            
        end
    end
    methods(Static) %math by hand
        % d matrix write by hand
        function SymExpr = d(j,m1,m2,seed)
            arguments
                j double{mustBePositive};
                m1 double
                m2 double
                seed sym
            end
            theta = seed;
            % m1 >= m2 set
            if (m1) < (m2)
                % d_{m^{\prime}, m}^{j}(-\beta)=d_{m, m^{\prime}}^{j}(\beta)=(-1)^{m^{\prime}-m} d_{m^{\prime}, m}^{j}(\beta)
                SymExpr = (-1)^(m2-m1)*Y_l__m.d(j,m2,m1,seed);
                return;
            end
            if (m1 <=0 ||(m1>0 && abs(m2)> m1))&& m2 <0 
                %d_{m^{\prime}, m}^{j}=(-1)^{m-m^{\prime}} d_{m, m^{\prime}}^{j}=d_{-m,-m^{\prime}}^{j}
                SymExpr = Y_l__m.d(j,-m2,-m1,seed);
                return;
            end
            % use property
            if isinteger(j) 
                if m1==0 && m2 == 0
                    % D_{0,0}^{\ell}(\alpha, \beta, \gamma)=d_{0,0}^{\ell}(\beta)=P_{\ell}(\cos \beta)
                    SymExpr = legendreP(j,cos(theta));
                    return
                end
                if m1 < 0 && m2 == 0
                    % d_{m 0}^{\ell}(\beta)=\sqrt{\frac{(\ell-m) !}{(\ell+m) !}} P_{\ell}^{m}(\cos \beta)
                    SymExpr = Y_l__m(j,m1,sqrt(4*sym(pi)/(2*j+1)),1).explicitformula(...
                        'vpa',false,'seed',[char(seed),'0']);
                    return;
                end
            end
            switch j
                case 1/2
                    if m1 ==1/2 && m2 ==1/2
                        SymExpr = cos(theta/2);
                    elseif m1 ==1/2 && m2 ==-1/2
                        SymExpr = -sin(theta/2);
                    else
                        SymExpr = sym(0);
                    end
                case 1
                    if m1 ==1 && m2 ==1
                        SymExpr = 0.5*(1+cos(theta));
                    elseif m1 ==1 && m2 ==-1
                        SymExpr = 0.5*(1-cos(theta));
                    elseif m1 ==1 && m2 ==0
                        SymExpr = -1/2^0.5*sin(theta);
                    elseif m1 ==0 && m2 ==0
                        SymExpr = cos(theta);
                    else
                        SymExpr = sym(0);
                    end
                case 3/2
                    switch m1
                        case 3/2
                            switch m2
                                case 3/2
                                    SymExpr = 1/2*(cos(theta)+1)*cos(theta/2);
                                case 1/2
                                    SymExpr = -sqrt(3)/2*(cos(theta)+1)*sin(theta/2);
                                case -1/2
                                    SymExpr = sqrt(3)/2*(-cos(theta)+1)*cos(theta/2);
                                case -3/2
                                    SymExpr = -1/2*(-cos(theta)+1)*sin(theta/2);
                            end
                        case 1/2
                            switch m2
                                case 1/2
                                    SymExpr = 1/2*(3*cos(theta)-1)*cos(theta/2);
                                case -1/2
                                    SymExpr = -1/2*(3*cos(theta)+1)*sin(theta/2);
                            end
                    end
                case 2
                    switch m1
                        case 2
                            switch m2
                                case 2
                                    SymExpr = 1/4*(cos(theta)+1)^2;
                                case 1
                                    SymExpr = -1/2*sin(theta)*(cos(theta)+1);
                                case -1
                                    SymExpr = sqrt(3)/2*(-cos(theta)+1)*cos(theta/2);
                                case -2
                                    SymExpr = -1/2*sin(theta)*(-cos(theta)+1);
                                case 0
                                    SymExpr = sqrt(3/8)*sin(theta)^2;
                            end
                        case 1
                            switch m2
                                case 1
                                    SymExpr = 1/2*(cos(2*theta)+cos(theta));
                                case -1
                                    SymExpr = 1/2*(-cos(2*theta)+cos(theta));
                                case 0
                                    SymExpr = -sqrt(3/8)*sin(2*theta);
                            end
                        case 0
                            switch m2
                                case 0
                                    SymExpr = 1/2*(3*cos(theta)^2-1);
                            end
                    end
                %case 5/2   
                %case 3    
                otherwise
                    SymExpr = Y_l__m.d_mm__j(j,m1,m2,seed);
            end
        end
        % Use jacobiP to define d metric(Wigner Origin format)
        function SymExpr = d_mm__j(j,m1,m2,seed)
            arguments
                j double{mustBePositive};
                m1 double
                m2 double
                seed sym
            end
            if m1 <0 &&  m2 <0
                %d_{m^{\prime}, m}^{j}=(-1)^{m-m^{\prime}} d_{m, m^{\prime}}^{j}=d_{-m,-m^{\prime}}^{j}
                SymExpr = Y_l__m.d_mm__j(j,-m2,-m1,seed);
                return;
            end
            theta = seed;
            % d_{m^{\prime} m}^{j}(\beta)=(-1)^{\lambda}\left(\begin{array}{c}
            % 2 j-k \\
            % k+a
            % \end{array}\right)^{\frac{1}{2}}\left(\begin{array}{c}
            % k+b \\
            % b
            % \end{array}\right)^{-\frac{1}{2}}\left(\sin \frac{\beta}{2}\right)^{a}\left(\cos \frac{\beta}{2}\right)^{b} P_{k}^{(a, b)}(\cos \beta)
            
            % https://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_(small)_d-matrix
            k = min([j+m2,j-m2,j+m1,j-m1]);
            switch k
                case j+m2
                    a = m1 - m2;lambda = m1-m2;
                case j-m2
                    a = m2 - m1;lambda = 0;
                case j+m1
                    a = m2 - m1;lambda = 0;
                case j-m1
                    a = m1 - m2;lambda = m1-m2;
            end
            b = 2*j-2*k-a;
            % nchoosek(n,k)
            % J = jacobiP(1,a,b,x)
            SymExpr = (-1)^lambda*...
                nchoosek(2*j-k,k+a)^0.5*nchoosek(k+b,b)^(-0.5)*...
                sin(theta/2)^a* cos(theta/2)^b*...
                jacobiP(k,a,b,cos(theta));
        end
        % The associated Legendre polynomials Plm can be defined recursively as:
        function PlmExpr = Plm(l,m,seed,options)
            % \begin{equation}
            % \begin{aligned}
            % P_{0}^{0}(\cos \theta) &=1 \\
            % P_{m}^{m}(\cos \theta) &=(-2l+1) \sin \theta P_{m-1}^{m-1}(\cos \theta) \\
            % P_{l=m+1}^{m}(\cos \theta) &=(2 l-1) \cos \theta P_{m}^{m}(\cos \theta) \\
            % P_{l}^{m}(\cos \theta) &=\frac{2 l-1}{l-m} \cos \theta P_{l-1}^{m}(\cos \theta)+\frac{1-l-m}{l-m} P_{l-2}^{m}(\cos \theta)
            % \end{aligned}
            % \end{equation}
            %
            arguments
                l double{mustBeNonnegative,mustBeInteger};
                m double{mustBeInteger}
                seed sym{} = (sym('theta','real'));
                options.ClosedForm = false;
                options.triangle = true;
            end
            optionsCell = namedargs2cell(options);
            if abs(m)>l
                ME = MException('vasplib:Y_l__m:WrongInput', ...
                    'Variable m:%d out of the range of [-l,l](l:%d)',m,l);
                throw(ME);
            end
            if strcmp(string(seed),"theta") && ~options.triangle
                x = sym('x','real');
                y = (1-x^2)^0.5;
            elseif options.triangle
                x  = cos(seed);
                y  = sin(seed);
            else
                x = seed;
                y = (1-x^2)^0.5; 
            end
            if options.ClosedForm
                %
                % \begin{equation}
                % P_{l}^{m}(x)=(-1)^{m} \cdot 2^{l} \cdot\left(1-x^{2}\right)^{m / 2} \cdot \sum_{k=m}^{l} \frac{k !}{(k-m) !} \cdot x^{k-m} \cdot\left(\begin{array}{l}
                % l \\
                % k
                % \end{array}\right)\left(\begin{array}{l}
                % \frac{l+k-1}{2} \\
                % l
                % \end{array}\right)
                % \end{equation}
                if  options.triangle
                    PlmExpr_pre = (-1)^m * 2^l * (y)^m;
                else
                    PlmExpr_pre = (-1)^m * 2^l * ((1-x^2)^(1/2))^m;
                end
                % binomialfunc = @Y_l__m.binomial
                binomialfunc = @nchoosek;
                k = m;
                PlmExpr_tail = factorial(k)/factorial(k-m)*x^(k-m)*binomialfunc(sym(l),k)*binomialfunc((l+k-1)/2,sym(l));
                for k = m+1:l
                PlmExpr_tail =PlmExpr_tail + factorial(k)/factorial(k-m)*x^(k-m)*binomialfunc(sym(l),k)*binomialfunc((l+k-1)/2,sym(l));
                end
                PlmExpr = simplify(PlmExpr_pre*simplify(PlmExpr_tail));
                return;
            end
            % P_{\ell}^{-m}=(-1)^{m} \frac{(\ell-m) !}{(\ell+m) !} P_{\ell}^{m}
            if m < 0
                PlmExpr = (-1)^(-m) * factorial(l+m)/factorial(l-m) * Y_l__m.Plm(-m,l,seed,optionsCell{:}); 
                return;
            end
            if l==0 && m ==0
                PlmExpr = 1;
                return;
            elseif l == m
                % first rule
                %PlmExpr = -(2*l-1)*y*Y_l__m.Plm(l-1,m-1,seed,optionsCell{:}); 
                PlmExpr = (-1)^l*Y_l__m.DFactorial(2*l-1)*y^l;
                return;
            elseif m == l-1
                % second rule
                PlmExpr = (2*m+1)*x*Y_l__m.Plm(m,m,seed,optionsCell{:});  
            else
                PlmExpr = (2*l-1)/(l-m) *x * Y_l__m.Plm(l-1,m,seed,optionsCell{:})...
                    + (1-l-m)/(l-m) * Y_l__m.Plm(l-2,m,seed,optionsCell{:});
            end
        end
        % fully normalized associated Legendre polynomials Nlm
        function NlmExpr = Nlm(l,m,seed,options)
            arguments
                l {mustBeNonnegative,mustBeInteger};
                m {mustBeInteger}
                seed sym{} = (sym('theta','real'));
                options.ClosedForm = false;
                options.triangle = true;
            end
            optionsCell = namedargs2cell(options);
            NlmExpr = (-1)^m * sqrt(((l+1/2))*factorial((l-m))/factorial((l+m)))* ...
                Y_l__m.Plm(l,m,seed,optionsCell{:});
        end
        % Ylm
        function YlmExpr = InnerY_l__m(l,m,seed1,seed2,options)
            arguments
                l double{mustBeNonnegative,mustBeInteger};
                m double{mustBeInteger}
                seed1 sym{} = (sym('theta','real'));
                seed2 sym{} = (sym('phi','real'));
                options.ClosedForm = false;
                options.triangle = true;
            end
            optionsCell = namedargs2cell(options);
            if options.triangle 
                phi  = seed2;
                YlmExpr = (-1)^m/sqrt(2*sym(pi))*...
                    Y_l__m.Nlm(l,m,seed1,optionsCell{:})*exp(1i*m*phi);
                return;
            else
                if strcmp(string(seed2),"phi") &&  strcmp(string(seed1),"theta")
                    x = sym('x','real');
                    y = sym('y','real');
                    z = sym('z','real');
                    r = sym('r','real');
                else
                    x = seed1(1); y = seed1(2); z =seed1(3);
                    r = seed2; 
                end
                YlmExpr = sqrt((2*l+1)*factorial(l-m)/factorial(l+m)/(4*sym(pi)))* ...
                    Y_l__m.Plm(l,m,z/r,optionsCell{:})*(1-(z/r)^2)^(-m/2)*(x/r+1i*y/r)^m;
                %YlmExpr = (-1)^m/sqrt(2*sym(pi))*...
                %    Y_l__m.Nlm(l,m,z/r,optionsCell{:})*(1-(z/r)^2)^(-m/2)*(x+1i*y)^m;
                YlmExpr = simplify(YlmExpr);
            end
        end
    end
    methods(Static) %math outer
        function DFact = DFactorial(n)
            % Double factorial function, n!!
            %  https://en.wikipedia.org/wiki/Double_factorial
            %
            % n!! = 1 for both n == 1 and n == 0.
            if isempty(n) || (numel(n) > 1) || (n < 0) || (mod(n,1) ~= 0)
                error('The sky is falling. n must be scalar, non-negative, integer.')
            end

            DFact = 1;
            if n > 1
                start = 1 + mod(n + 1,2); % caters for either parity of n
                DFact = prod(start:2:n);
            end
        end
        function c = binomial(n, k)
            % BINOMIAL Binomial coefficient.
            %
            %  If the arguments are both non-negative integers with 0 <= K <= N, then
            %  BINOMIAL(N, K) =  N!/K!/(N-K)!, which is the number of distinct sets of
            %  K objects that can be chosen from N distinct objects.
            %  When N or K(or both) are N-D matrices, BINOMIAL(N, K) is the coefficient
            %  for each pair of elements.
            %
            %  If N and K are integers that do not satisfy 0 <= K <= N, or N and K are
            %  non-integers , then the general definition is used, that is
            %
            %       BINOMIAL(N, K) = GAMMA(N+1) / (GAMMA(K+1) / GAMMA(N-K+1))
            %
            %  If N is a non-negative integer, BINOMIAL(N) = BINOMIAL(N, 0:N)
            %
            %  BINOMIAL only supports floating point input arguments.
            % Mukhtar Ullah
            % mukhtar.ullah@informatik.uni-rostock.de
            % Created: 5 Oct, 2004
            % Last Modified: 18 Oct, 2010
            error(nargchk(1, 2, nargin));
            if nargin == 1
                validateattributes(n, {'numeric'}, {'scalar', 'real', 'nonnegative'});
                c = diag(fliplr(pascal(floor(n) + 1))).';
            else
                assert(isscalar(n) || isscalar(k) || isequal(size(n),size(k)), ...
                    'Non-scalar arguments must have the same size.')
                validateattributes([n(:); k(:)], {'numeric'}, {'real', 'nonnegative'});
                c = exp(gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1));  % binomial coefficient
                i = ( n==floor(n+.5) & k==floor(k+.5) );
                c(i) = floor(c(i)+.5);                                  % number of combinations
            end
        end
        function cg = CG(j1,m1,j2,m2,j,m)
            assert(isscalar(j1) && isscalar(m1) && isscalar(j2) && isscalar(m2) && isscalar(j) && isscalar(m),'All inputs must be scalars.')
            %---check conditions---%
            if j1 < 0 || j2 < 0 || j < 0 || ...
                    mod(2*j1,1) ~= 0 || mod(2*j2,1) ~= 0 || mod(2*j,1) ~= 0 || ...
                    mod(2*m1,1) ~= 0 || mod(2*m2,1) ~= 0 || mod(2*m,1) ~= 0 || ...
                    abs(m1) > j1 || abs(m2) > j2 || abs(m) > j || ...
                    j1+m1 < 0 || j2+m2 < 0 || j+m < 0 || j1+j2+j < 0 ||...
                    mod(j1+m1,1) ~= 0 || mod(j2+m2,1) ~= 0 || mod(j+m,1) ~= 0 || ...
                    mod(j1+j2+j,1) ~= 0
                
                error(sprintf('Clebsch-Gordan coefficient only defined if: \n 1. j1, j2, j are integer or half-integer non-negative numbers. \n 2. m1, m2, m are integer or half-integer numbers. \n 3. abs(m1)<=j1, abs(m2)<=j2, abs(m)<=j \n 4. j1+m1, j2+m2, j+m, j1+j2+j are integer non-negative numbers.')) %#ok<SPERR>
                
            elseif m1+m2-m ~= 0 || j < abs(j1-j2) || j > j1+j2
                
                C = 0;
                return
                
            end
            
            %---compute valid k values for summation---%
            k = max([0,j2-j-m1,j1-j+m2]):min([j1+j2-j,j1-m1,j2+m2]);
            %---check for stability---%
            if j+j1-j2 > 21 || j+j2-j1 > 21 || j1+j2-j > 21 || j1+j2+j+1 > 21 || ...
                    j+m > 21 || j-m > 21 || j1+m1 > 21 || j1-m1 > 21 || ...
                    j2+m2 > 21 || j2-m2 > 21 || any(k > 21) || ...
                    any(j1+j2-j-k > 21) || any(j1-m1-k > 21) || ...
                    any(j2+m2-k > 21) || any(j-j2+m1+k > 21) || any(j-j1-m2+k > 21)
                
                warning('The argument to one or more of the factorials used in the computation of the requested Clebsch-Gordan coefficient is greater than 21, this can result in inaccuracies (see Matlab documentation for the FACTORIAL function).') %#ok<WNTAG>
                
            end
            %---compute coefficient---%
            cg = sqrt((2*j+1)*factorial(j+j1-j2)*factorial(j+j2-j1)*factorial(j1+j2-j)/factorial(j1+j2+j+1))*...
                sqrt(factorial(j+m)*factorial(j-m)*factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2))*...
                sum(((-1).^k)./(factorial(k).*factorial(j1+j2-j-k).*factorial(j1-m1-k).*factorial(j2+m2-k).*factorial(j-j2+m1+k).*factorial(j-j1-m2+k)));
        end
        function  W = w3j(j1,j2,j3,m1,m2,m3)
            %-------------------------------------------------------------------------%
            %Filename:  w3j.m
            %Author:    Oliver Johnson
            %Date:      6/7/2011
            %
            % W3J computes the Wigner 3-j coefficients.
            %
            % Inputs:
            %   j1 - A scalar giving the first total angular momentum.
            %   m1 - A scalar giving the projection of the first total angular
            %        momentum.
            %   j2 - A scalar giving the second total angular momentum.
            %   m2 - A scalar giving the projection of the second total angular
            %        momentum.
            %   j3 - A scalar giving the third total angular momentum.
            %   m3 - A scalar giving the projection of the third total angular
            %        momentum.
            %
            % Outputs:
            %   W - A scalar giving the required Wigner 3-j coefficient.
            %
            % [1] http://en.wikipedia.org/wiki/3-jm_symbol.
            %-------------------------------------------------------------------------%
            W = (-1)^(j1-j2-m3)*Y_l__m.CG(j1,m1,j2,m2,j3,-m3)/sqrt(2*j3+1);
        end
        function  W = w6j(a,b,c,d,e,f)
            %-------------------------------------------------------------------------%
            %Filename:  w6j.m
            %Author:    Oliver Johnson
            %Date:      6/7/2011
            %
            % W6J computes the Wigner 6-j coefficients.
            %
            % Inputs:
            %   a,b,c,d,e,f - Scalar arguments to the 6-j symbol from the form:
            %
            %                                   {a b c}
            %                                   {d e f}
            %
            % Outputs:
            %   W - A scalar giving the required Wigner 6-j coefficient.
            %
            % Note: Tested on all 2264 cases given at [3], and results were within eps.
            %
            % [1] Varshoalovich; Quantum Theory of Angular Momentum. (1988). p. 293.
            % [2] Messiah; Quantum Mechanics, Vol. 2. (1962). p.1063.
            % [3] http://www.strw.leidenuniv.nl/~mathar/progs/6jSymb
            % [4] Format of nlims taken from:
            %     http://www.mathworks.com/matlabcentral/fileexchange/20619
            %-------------------------------------------------------------------------%
            assert(isscalar(a) && isscalar(b) && isscalar(c) && isscalar(d) && isscalar(e) && isscalar(f),'All inputs must be scalars.')
            %---check conditions---%
            if ~Y_l__m.triangular_cond(a,b,c) || ~Y_l__m.triangular_cond(c,d,e) || ...
                    ~Y_l__m.triangular_cond(a,e,f) || ~Y_l__m.triangular_cond(b,d,f) || ...
                    mod(a+b+c,1) ~= 0 || mod(c+d+e,1) ~= 0 || ...
                    mod(a+e+f,1) ~= 0 || mod(b+d+f,1) ~= 0
                
                W = 0;
                return
                
            end
            
            %---compute valid n values for summation---%
            nlims(1) = 0;
            nlims(2) = a+b+c;
            nlims(3) = c+d+e;
            nlims(4) = a+e+f;
            nlims(5) = b+d+f;
            nlims(6) = a+b+d+e;
            nlims(7) = a+c+d+f;
            nlims(8) = b+c+e+f;
            n = max(nlims(1:5)):min(nlims(6:8));
            
            %---check for stability---%
            if a+b-c > 21 || a-b+c > 21 || -a+b+c > 21 || a+b+c+1 > 21 || ...
                    c+d-e > 21 || c-d+e > 21 || -c+d+e > 21 || c+d+e+1 > 21 || ...
                    a+e-f > 21 || a-e+f > 21 || -a+e+f > 21 || a+e+f+1 > 21 || ...
                    b+d-f > 21 || b-d+f > 21 || -b+d+f > 21 || b+d+f+1 > 21 || ...
                    any(n+1 > 21) || ...
                    any(n-nlims(2) > 21) || any(n-nlims(3) > 21) || ...
                    any(n-nlims(4) > 21) || any(n-nlims(5) > 21) || ...
                    any(nlims(6)-n > 21) || any(nlims(7)-n > 21) || ...
                    any(nlims(8)-n > 21)
                
                warning('The argument to one or more of the factorials used in the computation of the requested Clebsch-Gordan coefficient is greater than 21, this can result in inaccuracies (see Matlab documentation for the FACTORIAL function).') %#ok<WNTAG>
                
            end
            
            %---compute coefficient---%
            W = Y_l__m.del(a,b,c)*del(c,d,e)*Y_l__m.del(a,e,f)*del(b,d,f)*...
                sum((((-1).^n).*factorial(n+1))./...
                (factorial(n-nlims(2)).*factorial(n-nlims(3)).*...
                factorial(n-nlims(4)).*factorial(n-nlims(5)).*...
                factorial(nlims(6)-n).*factorial(nlims(7)-n).*...
                factorial(nlims(8)-n)));
            
        end
        function  W = w9j(a,b,c,d,e,f,g,h,j)
            %-------------------------------------------------------------------------%
            %Filename:  w9j.m
            %Author:    Oliver Johnson
            %Date:      6/7/2011
            %
            % W9J computes the Wigner 9-j coefficients.
            %
            % Inputs:
            %   a,b,c,d,e,f,g,h,j - Scalar arguments to the 9-j symbol from the form:
            %
            %                                   {a b c}
            %                                   {d e f}
            %                                   {g h j}
            %
            % Outputs:
            %   W - A scalar giving the required Wigner 9-j coefficient.
            %
            % Note: Tested on first 87 cases from Table 10.13 of Ref. [1].
            %
            % [1] Varshoalovich; Quantum Theory of Angular Momentum. (1988). p. 340.
            % [2] Format of xlims taken from:
            %     http://www.mathworks.com/matlabcentral/fileexchange/20619
            % [3] Weisstein, Eric W. "Triangular Inequalities."
            %     From MathWorld--A Wolfram Web Resource.
            %     http://mathworld.wolfram.com/TriangularInequalities.html
            %-------------------------------------------------------------------------%
            assert(isscalar(a) && isscalar(b) && isscalar(c) && ...
                isscalar(d) && isscalar(e) && isscalar(f) && ...
                isscalar(g) && isscalar(h) && isscalar(j),'All inputs must be scalars.')
            
            %---check conditions---%
            args = [a b c d e f g h j];
            if any(args < 0) || any(mod(2*args,1) ~= 0)
                
                error('All arguments to Wigner 9-j symbol must be integer or half-integer non-negative numbers.')
                
            elseif ~triangular_cond(a,b,c) || ~triangular_cond(d,e,f) || ...
                    ~triangular_cond(g,h,j) || ~triangular_cond(a,d,g) || ...
                    ~triangular_cond(b,e,h) || ~triangular_cond(c,f,j)
                
                W = 0;
                return
                
            end
            
            %---compute valid x values for summation---%
            xlims(1) = a+j;
            xlims(2) = b+f;
            xlims(3) = d+h;
            xlims(4) = a-j;
            xlims(5) = b-f;
            xlims(6) = d-h;
            
            %---compute W---%
            W = 0;
            for x = max(xlims(4:6)):min(xlims(1:3))
                
                W = W + (-1)^(2*x)*(2*x+1)*...
                    w6j(a,b,c,f,j,x)*w6j(d,e,f,b,x,h)*w6j(g,h,j,x,a,d);
                
            end
            
        end
        function cg = ClebschGordan(j1,j2,j,m1,m2,m)
            % ClebschGordan.m by David Terr, Raytheon, 6-17-04
            % Modified on 11-9-04
            
            % ClebschGordan(j1,j2,j,m1,m2,m) returns the Clebsch-Gordan coefficient <j1,j2,m1,m2|j1,j2,j,m>.
            % This program requires first downloading Wigner3j.m.
            % error checking
            if ( 2*j1 ~= floor(2*j1) || 2*j2 ~= floor(2*j2) || 2*j ~= floor(2*j) ...
                    || 2*m1 ~= floor(2*m1) || 2*m2 ~= floor(2*m2) || 2*m ~= floor(2*m) )
                error('All arguments must be integers or half-integers.');
                return;
            end
            
            if m1 + m2 ~= m
                warning('m1 + m2 must equal m.');
                cg = 0;
                return;
            end
            
            if ( j1 - m1 ~= floor ( j1 - m1 ) )
                warning('2*j1 and 2*m1 must have the same parity');
                cg = 0;
                return;
            end
            
            if ( j2 - m2 ~= floor ( j2 - m2 ) )
                warning('2*j2 and 2*m2 must have the same parity');
                cg = 0;
                return;
            end
            
            if ( j - m ~= floor ( j - m ) )
                warning('2*j and 2*m must have the same parity');
                cg = 0;
                return;
            end
            
            if j > j1 + j2 || j < abs(j1 - j2)
                warning('j is out of bounds.');
                cg = 0;
                return;
            end
            
            if abs(m1) > j1
                warning('m1 is out of bounds.');
                cg = 0;
                return;
            end
            
            if abs(m2) > j2
                warning('m2 is out of bounds.');
                cg = 0;
                return;
            end
            
            if abs(m) > j
                warning('m is out of bounds.');
                cg = 0;
                return;
            end
            
            cg = (-1)^(j1-j2+m) * sqrt(2*j + 1) * Y_l__m.Wigner3j(j1,j2,j,m1,m2,-m);
            
            
            % Reference: Clebsch-Gordan Coefficient entry of Eric Weinstein's Mathworld: http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html
        end
        function  W = Wigner3j( j123, m123 )
            % Compute the Wigner 3j symbol using the Racah formula.
            %
            % W = Wigner3j( J123, M123 )
            %
            % J123 = [J1, J2, J3].
            % M123 = [M1, M2, M3].
            % All Ji's and Mi's have to be integeres or half integers (correspondingly).
            %
            % According to seletion rules, W = 0 unless:
            %   |Ji - Jj| <= Jk <= (Ji + Jj)    (i,j,k are permutations of 1,2,3)
            %   |Mi| <= Ji    (i = 1,2,3)
            %    M1 + M2 + M3 = 0
            %
            % Reference:
            % Wigner 3j-Symbol entry of Eric Weinstein's Mathworld:
            % http://mathworld.wolfram.com/Wigner3j-Symbol.html
            %
            % Inspired by Wigner3j.m by David Terr, Raytheon, 6-17-04
            %  (available at www.mathworks.com/matlabcentral/fileexchange).
            %
            % By Kobi Kraus, Technion, 25-6-08.
            % Updated 1-8-13.
            
            j1 = j123(1); j2 = j123(2); j3 = j123(3);
            m1 = m123(1); m2 = m123(2); m3 = m123(3);
            
            % Input error checking
            if any( j123 < 0 ),
                error( 'The j must be non-negative' )
            elseif any( rem( [j123, m123], 0.5 ) ),
                error( 'All arguments must be integers or half-integers' )
            elseif any( rem( (j123 - m123), 1 ) )
                error( 'j123 and m123 do not match' );
            end
            
            % Selection rules
            if ( j3 > (j1 + j2) ) || ( j3 < abs(j1 - j2) ) ... % j3 out of interval
                    || ( m1 + m2 + m3 ~= 0 ) ... % non-conserving angular momentum
                    || any( abs( m123 ) > j123 ), % m is larger than j
                W = 0;
                return
            end
            
            % Simple common case
            if ~any( m123 ) && rem( sum( j123 ), 2 ), % m1 = m2 = m3 = 0 & j1 + j2 + j3 is odd
                W = 0;
                return
            end
            
            % Evaluation
            t1 = j2 - m1 - j3;
            t2 = j1 + m2 - j3;
            t3 = j1 + j2 - j3;
            t4 = j1 - m1;
            t5 = j2 + m2;
            
            tmin = max( 0,  max( t1, t2 ) );
            tmax = min( t3, min( t4, t5 ) );
            
            t = tmin : tmax;
            W = sum( (-1).^t .* exp( -ones(1,6) * gammaln( [t; t-t1; t-t2; t3-t; t4-t; t5-t] +1 ) + ...
                gammaln( [j1+j2+j3+1, j1+j2-j3, j1-j2+j3, -j1+j2+j3, j1+m1, j1-m1, j2+m2, j2-m2, j3+m3, j3-m3] +1 ) ...
                * [-1; ones(9,1)] * 0.5 ) ) * (-1)^( j1-j2-m3 );
            
            % Warnings
            if isnan( W )
                warning( 'MATLAB:Wigner3j:NaN', 'Wigner3J is NaN!' )
            elseif isinf( W )
                warning( 'MATLAB:Wigner3j:Inf', 'Wigner3J is Inf!' )
            end
        end
        function  W = Wigner3j_sym(j123, m123)
                % Computes exact value of Wigner 3-j symbol with symbolic toolbox
                % Tested principally for integral input parameters
                % Does not check for ~(integer || half-integer)
                % Includes additional symbolic function Y_l__m.sym_fact() for exact factorials
                % Note that applying eval() to function value gives floating-point approximation
                a = j123(1); b = j123(2); c = j123(3);
                alf = m123(1); bet = m123(2); gam = m123(3);
                % apply selection rules
                tj_iszero = false;
                if abs(alf) > abs(a) || abs(bet) > abs(b) || abs(gam) > abs(c)
                    tj_iszero = true;
                end
                if alf + bet + gam ~= 0
                    tj_iszero = true;
                end
                if abs(c) > abs(a + b) || abs(c) < abs(a - b)
                    tj_iszero = true;
                end
                
                % use Racah formula
                if ~tj_iszero
                    % compute pre-factors
                    tri_coef = Y_l__m.sym_fact(a + b - c)*Y_l__m.sym_fact(a - b + c)*Y_l__m.sym_fact(-a + b + c)/Y_l__m.sym_fact(a + b + c + 1);
                    pre_coef = Y_l__m.sym_fact(a + alf)*Y_l__m.sym_fact(a - alf)...
                        *Y_l__m.sym_fact(b + bet)*Y_l__m.sym_fact(b - bet)...
                        *Y_l__m.sym_fact(c + gam)*Y_l__m.sym_fact(c - gam);
                    
                    % find maximum possible number of terms (assuming t increases from zero), then add up series
                    tmax = max([a+alf,a-alf, b+bet,b-bet, c+gam,c-gam,...
                        a+b-c,b+c-a,c+a-b]);
                    t_sum = sym(0);
                    for t = 0:tmax
                        
                        % apply termwise selection rule, sum over allowed terms
                        if c-b+t+alf >= 0 && c-a+t-bet >=0 && a+b-c-t >= 0 && a-t-alf >=0 && b-t+bet >= 0
                            t_denom = Y_l__m.sym_fact(t)*Y_l__m.sym_fact(a+b-c-t)*...
                                Y_l__m.sym_fact(c-b+t+alf)*Y_l__m.sym_fact(c-a+t-bet)*...
                                Y_l__m.sym_fact(a-t-alf)*Y_l__m.sym_fact(b-t+bet);
                            t_sum = t_sum + sym((-1)^t)/t_denom;
                        end
                    end
                    
                    W = sym((-1)^(a-b-gam))*sqrt(tri_coef*pre_coef)*t_sum;
                end
                
                % if zero
                if tj_iszero
                    W = sym(0);
                end
                % use simple() to yield clean result
                W = simplify(W);
        end
    end
    %% tool
    methods(Static) % addtional tool
        function OutExpr = nlm2atomic(l,m,n,options)
            arguments
                l = 1;
                m = 0; % m for real spherical harmonic see <https://en.wikipedia.org/wiki/Table_of_spherical_harmonics>
                n = 0;
                options.outputformat = 'sym';
            end
            L_dist = containers.Map([0,1,2,3,4,5,6],["s","p","d","f","g","h","i"]);
            if strcmp(options.outputformat ,'sym')
                if n ==0
                    str1 = "";
                    str2 = "";
                else
                    str1 = num2str(n);
                    str2 = num2str(n);
                end
                if l >= 0 && l <=6
                    str1 = str1 + L_dist(l)+"_" + Y_l__m.lm2str(l,abs(m)); % m > 0 component
                    str2 = str2 + L_dist(l)+"_" + Y_l__m.lm2str(l,-abs(m));
                elseif l < 0

                else

                end
                if m == 0
                    OutExpr = sym(str1,'real');
                else
                    % See: https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
                    % spherical harmonocs expanded by tesseral spherical
                    % harmonics
                    % $$
                    % Y_{\ell}^{m}= \begin{cases}\frac{1}{\sqrt{2}}\left(Y_{\ell|m|}-i Y_{\ell,-|m|}\right) & \text { if } m<0 \\ Y_{\ell 0} & \text { if } m=0 \\ \frac{(-1)^{m}}{\sqrt{2}}\left(Y_{\ell|m|}+i Y_{\ell,-|m|}\right) & \text { if } m>0\end{cases}
                    % $$
                    % OutExpr = -1*sign(m)*sym(str1,'real')/sqrt(2) + (-1)^m*1i/sqrt(2)*sym(str2,'real');
                    if m < 0
                        OutExpr = 1/sqrt(2) * (sym(str1,'real')-1i*sym(str2,'real'));
                    else
                        OutExpr = (-1)^m/sqrt(2) * (sym(str1,'real')+1i*sym(str2,'real'));
                    end
                end
            end
        end
        function Outstr = l__m2str(l,m)
            if l <1

            end
            switch l
                case 0
                    switch m
                        case 0
                            Outstr = "s";
                    end
                case 1
                    switch m
                        case 0
                            Outstr = "z";
                        case -1
                            Outstr = "y";
                        case 1
                            Outstr = "x";
                    end
                case 2
                    switch m
                        case 0
                            Outstr = "z2";
                        case -1
                            Outstr = "yz";
                        case 1
                            Outstr = "xz";
                        case -2
                            Outstr = "xy";
                        case 2
                            Outstr = "x2my2";
                    end
                case 3
                    switch m
                        case 0
                            Outstr = "z3";
                        case -1
                            Outstr = "yz2";
                        case 1
                            Outstr = "xz2";
                        case -2
                            Outstr = "xyz";
                        case 2
                            Outstr = "zx2my2";
                        case -3
                            Outstr = "3x2ymy3";
                        case 3
                            Outstr = "x2m3y2";
                    end
                otherwise
                    Outstr = num2str(abs(m));
                    if sign(m) == -1
                        Outstr = [Outstr,'__bar'];
                    end
            end


        end
        function Outstr = lm2str(l,m)
            if l <1

            end
            switch l
                case 0
                    switch m
                        case 0
                            Outstr = "";
                    end
                case 1
                    switch m
                        case 0
                            Outstr = "z";
                        case -1
                            Outstr = "y";
                        case 1
                            Outstr = "x";
                    end
                case 2
                    switch m
                        case 0
                            Outstr = "z2";
                        case -1
                            Outstr = "yz";
                        case 1
                            Outstr = "xz";
                        case -2
                            Outstr = "xy";
                        case 2
                            Outstr = "x2my2";
                    end
                case 3
                    switch m
                        case 0
                            Outstr = "z3";
                        case -1
                            Outstr = "yz2";
                        case 1
                            Outstr = "xz2";
                        case -2
                            Outstr = "xyz";
                        case 2
                            Outstr = "zx2my2";
                        case -3
                            Outstr = "3x2ymy3";
                        case 3
                            Outstr = "x2m3y2";
                    end
                otherwise
                    Outstr = num2str(abs(m));
                    if sign(m) == -1
                        Outstr = [Outstr,'__bar'];
                    end
            end
        end
        % find duplicate entries in the list
        function [dupNames, dupNdxs] = getDuplicates(aList)
            [uniqueList,~,uniqueNdx] = unique(aList);
            N = histc(uniqueNdx,1:numel(uniqueList));
            dupNames = uniqueList(N>1);
            dupNdxs = arrayfun(@(x) find(uniqueNdx==x), find(N>1), ...
                'UniformOutput',false);
        end
        function Ind = IndexOfMultiples(A)
            T   = true(size(A));
            off = false;
            A   = A(:);
            for iA = 1:numel(A)
                if T(iA)          % if not switched already
                    d = (A(iA) == A);
                    if sum(d) > 1   % More than 1 occurrence found
                        T(d) = off;  % switch all occurrences
                    end
                end
            end
            Ind = find(~T);
        end
        function T = isMultiple(A)
            % T = isMultiple(A)
            % INPUT:  A: Numerical or CHAR array of any dimensions.
            % OUTPUT: T: TRUE if element occurs multiple times anywhere in the array.
            %
            % Tested: Matlab 2009a, 2015b(32/64), 2016b, 2018b, Win7/10
            % Author: Jan, Heidelberg, (C) 2021
            % License: CC BY-SA 3.0, see: creativecommons.org/licenses/by-sa/3.0/
            
            T        = false(size(A));
            [S, idx] = sort(A(:).');
            m        = [false, diff(S) == 0];
            if any(m)        % Any equal elements found:
                m(strfind(m, [false, true])) = true;
                T(idx) = m;   % Resort to original order
            end
        end
    end
    methods(Static)
        function tri = del(a,b,c)
            
            tri = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1)); %Eq. 8.2(1)
            
        end
        function tf = triangular_cond(a,b,c)
            
            tf = (c >= abs(a-b)) & (c <= a+b);
            
        end
        function sym_fact = sym_fact(n)
            % computes exact integer factorial
            if n == 0
                sym_fact = sym(1);
            end
            if n > 0
                sym_fact = prod(sym(1):sym(n));
            end
            if n < 0
                sym_fact = sym(inf);
            end
        end
    end
end
