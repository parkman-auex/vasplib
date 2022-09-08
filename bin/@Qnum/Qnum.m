classdef Qnum < HollowKnight
    %Qnums  Quantum numbers of a state
%   See also LEGENDRE.
    properties %(Access = private)
        q = 0;
        AN = 6;    
        n = 2;
        l = 1;
        m = 0;
        s = 0;
        sz = 0;
    end
    properties
        %coe =1;
    end
    properties(Dependent)
        hollow;
    end
    %% construction
    methods
        function Qnums = Qnum(AN,n,l,m,s,sz,q,coe,propArgs,options)
            arguments
                AN  = 6;
                n   = 2;%... Radius(n?)
                l   = 1;
                m   = 0;%... Magnetic quantum number
                s   = 0;
                sz  = 0;
                q   = 0;%... Azimuthal quantum number
                coe = 1;
                propArgs.?Qnum;
                options.set = false;
            end
            %
            if options.set
                q  = propArgs.q;
                AN  = propArgs.AN;
                n  = propArgs.n;
                l  = propArgs.l;
                m  = propArgs.m;
                s  = propArgs.s;
                sz = propArgs.sz;
                coe = propArgs.coe;
            else
                
            end
            %PropertyCell = namedargs2cell(options);
            Qnums.q = q;
            Qnums.AN = AN;
            Qnums.n = n;
            Qnums.l = l;
            Qnums.n = m;
            Qnums.s = s;
            Qnums.sz = sz;
            Qnums.coe = coe;


            % multi
            [nS_row,nS_col] = size(l);
            [nSz_row,nSz_col] = size(m);
            [ncoe_row,ncoe_col] = size(coe);
            [nn,nn_col] = size(n);
            nQnums = max([nS_row ,nSz_row, ncoe_row,nn]);
            nQnums_col = max([nS_col ,nSz_col, ncoe_col nn_col] );
            if nQnums == 1 && nQnums_col == 1
                return;
            else
                Qnums = Qnum.QnumL(AN,n,l,m,s,sz,q,coe,'nrow',nQnums,'ncol',nQnums_col);
            end
        end
    end
    methods(Static)
        function [Qnums,sL,szL] = QnumL(ANL,nL,lL,mL,sL,szL,qL,coeL,options)
            arguments
                ANL  = 6;
                nL   = 2;%... Radius(n?)
                lL   = 1;
                mL   = 0;%... Magnetic quantum number
                sL   = 0;
                szL  = 0;
                qL   = 0;%... Azimuthal quantum number
                coeL = 1;
                options.list = 0;
                options.nQnums = [];
                options.nQnums_col = [];
            end
            optionsCell = namedargs2cell(options);
            if options.list
                Qnums = Qnum(ANL,nL,lL,mL,sL,szL,qL,coeL);
                Qnums = repmat(Qnums,[options.list  1]);
                return;
            end
            if isa(ANL,'vasplib')
                quantumL = ANL.quantumL;
                ANL = ANL.elementL;
                [Qnums,sL,szL] = Qnum.QnumL(ANL,quantumL,optionsCell{:});
                return;
            end
            if ~isvector(ANL) && isscalar(ANL) && isvector(nL)
                coeL = nL;
                quantumL = ANL;
                ANL = quantumL(:,1);
                nL = quantumL(:,2);
                lL = quantumL(:,3);
                mL = quantumL(:,4);
                sL = quantumL(:,5);
                szL = quantumL(:,6);
                qL = quantumL(:,7);
            end
            if size(nL,2)>1
                quantumL = nL;
                nL = quantumL(:,1);
                lL = quantumL(:,2);
                mL = quantumL(:,3);
                if size(quantumL,2)>4
                    sL = quantumL(:,4);
                    szL = quantumL(:,5);
                else
                    szL = quantumL(:,4);
                    if find(sign(szL) == -1)
                        % SOC
                        szL = 1/2*sign(szL);
                    else
                        szL = zeros(size(szL));
                    end
                    sL = abs(szL);
                end
                [Qnums,sL,szL] = Qnum.QnumL(ANL,nL,lL,mL,sL,szL,coeL,optionsCell{:});
                return;
            end
            if isempty(options.nQnums)
                nQnums = size(mL,1);
            else
                nQnums = options.nQnums;
            end
            if isempty(options.nQnums_col)
                nQnums_col = 1;
            else
                nQnums_col = options.nQnums_col;
            end
            Qnums = Qnum(ANL(1),nL(1),lL(1),mL(1),sL(1),szL(1),qL(1),coeL(1));
            Qnums = repmat(Qnums,[nQnums  nQnums_col]);
            [ANL,nL,lL,mL,sL,szL,qL,coeL] = Qnum.StanderdInput(ANL,nL,lL,mL,sL,szL,qL,coeL,'nrow',nQnums,'ncol',nQnums_col);
            for i = 1:numel(Qnums)
                Qnums(i).AN = ANL(i);
                Qnums(i).n  = nL(i);
                Qnums(i).l  = lL(i);
                Qnums(i).m  = mL(i);
                Qnums(i).s  = sL(i);
                Qnums(i).sz = szL(i);
                Qnums(i).q  = qL(i);
                Qnums(i).coe  = coeL(i);
            end
        end

        %function [ANL,nL,lL,mL,sL,szL,qL,] = StanderdInput(ANL,nL,lL,mL,sL,szL,qL,nQnums,nQnums_col)
    end
    methods % reload

    end
    %% rotation
    methods
        function QmumObj = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
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
            %QmumObj = A.HollowMe();
            %%%%%%%%%%% GLOBAL COE %%%%%%%%%%%
            Coe_global = A.coe;
            %%%%%%%%%%% Spin Rotation %%%%%%%%%%%
            szL = [];
            if isa(Coe_global,'sym')
                scoeL = sym([]);
            else
                scoeL = [];
            end
            As = A.s;
            WDpreAs1 = exp(1i*RightorLeft*A.sz*alpha);
            %             WDpreA2 = coeL_A(2)*exp(1i*RightorLeft*( A.m)*alpha);
            if As == 0
                szL = 0;
                scoeL = [scoeL;1];
            else
                for szB = -As:1:As
                    WignerD_single_element_sz = (Y_l__m.d(As,A.sz,szB,beta));
                    szL = [szL;szB];
                    scoeL = [scoeL;WDpreAs1*WignerD_single_element_sz*exp(1i*RightorLeft*szB*gamma)];
                    %if A.parity == 1 suppose parity = 1
                    %    Ai.coe = Ai.coe;
                    %else
                    %    Ai.coe = -1*Ai.coe;
                    %end
                    % for timereversal?
                end
                [szL,scoeL] = HollowKnight.generalcontractrow(szL,scoeL);
            end
            if conjugate
                % timereveersal?
                scoeL = (-1).^(As - szL).* scoeL;
                szL = -szL;
            end
            %%%%%%%%%%% Orbital Rotation %%%%%%%%%%% 
            Y_lmTmp = Y_lm(A.l,A.m); % only one |l,m>
            Y_lmTmp_rotate = Y_lmTmp(1).rotateinner(abc,RightorLeft,immproper,conjugate,antisymmetry);
            if length(Y_lmTmp)>1
                Y_lmTmp_rotate = [Y_lmTmp_rotate,Y_lmTmp(2).rotateinner(abc,RightorLeft,immproper,conjugate,antisymmetry)];
            end
            Tesseralexpansion = Tesseral(Y_lmTmp_rotate);
            Al = Tesseralexpansion{1,1}(:,1);
            mL = Tesseralexpansion{1,1}(:,2);
            coeL = Tesseralexpansion{1,2};
            %%%%%%%%%% complex to real SH %%%%%%%%%%%
            % bad coding remove later
            %
            %             mL = [];
            %             if isa(Coe_global,'sym')
            %                 coeL = sym([]);
            %             else
            %                 coeL = [];
            %             end
            %             Al = A.l;
            %             if Al == 0
            %                 mL = 0;
            %                 coeL = [coeL;1];
            %             else
            %                 % A % Y_lm = Y_l__(-m)  Y_l__(m)
            %                 if A.m < 0
            %                     coeL_A(1)  =  1i/sqrt(2);
            %                     coeL_A(2)  =  -1i/sqrt(2)*(-1)^A.m;
            %                 elseif A.m == 0
            %                     coeL_A(1)  =  1/2;
            %                     coeL_A(2)  =  1/2;
            %                 elseif A.m > 0
            %                     coeL_A(1)  =  1/sqrt(2);
            %                     coeL_A(2)  =  1/sqrt(2)*(-1)^A.m;
            %                 end
            %                 absm = abs(A.m);
            %                 WDpreA1 = coeL_A(1)*exp(1i*RightorLeft*(-absm)*alpha);
            %                 WDpreA2 = coeL_A(2)*exp(1i*RightorLeft*( absm)*alpha);
            %                 %
            %                 for mB = -Al:1:Al
            %                     %                     % back to complexSH
            %                     %                     % See: https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
            %                     %                     % spherical harmonocs expanded by tesseral spherical
            %                     %                     % harmonics
            %                     %                     % $$
            %                     %                     % Y_{\ell}^{m}= \begin{cases}\frac{1}{\sqrt{2}}\left(Y_{\ell|m|}-i Y_{\ell,-|m|}\right) & \text { if } m<0 \\ Y_{\ell 0} & \text { if } m=0 \\ \frac{(-1)^{m}}{\sqrt{2}}\left(Y_{\ell|m|}+i Y_{\ell,-|m|}\right) & \text { if } m>0\end{cases}
            %                     %                     % $$
            %                     %                     if mB < 0
            %                     %                         coeL_B(1,1)  =  1i/sqrt(2);% 1/sqrt(2);
            %                     %                         coeL_B(2,1)  =  -1i/sqrt(2)*(-1)^mB;%-1i/sqrt(2);
            %                     %                     elseif mB == 0
            %                     %                         coeL_B(1,1)  =  1/2;
            %                     %                         coeL_B(2,1)  =  1/2;
            %                     %                     elseif mB > 0
            %                     %                         coeL_B(1,1)  =  1/sqrt(2);%1/sqrt(2)*(-1)^mB;
            %                     %                         coeL_B(2,1)  =  1/sqrt(2)*(-1)^mB;%1i/sqrt(2)*(-1)^mB;
            %                     %                     end
            %                     %                     absmB = abs(mB);
            %                     %                     % optimize thesecodes after comfirmation!
            %                     %                     mL = [mL;-absmB;absmB;-absmB;absmB;]; % Y_lm = Y_l__(-m)  Y_l__(m)
            %                     %                     WignerD_single_element1 = (Y_l__m.d(Al,-absm,-mB,beta));
            %                     %                     WignerD_single_element2 = (Y_l__m.d(Al,-absm,mB,beta));
            %                     %                     WignerD_single_element3 = (Y_l__m.d(Al, absm,-mB,beta));
            %                     %                     WignerD_single_element4 = (Y_l__m.d(Al, absm,mB,beta));
            %                     %                     coeL_C(1,1) = WDpreA1*WignerD_single_element1*exp(1i*RightorLeft*mB*gamma)*coeL_B(1,1);
            %                     %                     coeL_C(2,1) = WDpreA1*WignerD_single_element2*exp(1i*RightorLeft*mB*gamma)*coeL_B(2,1);
            %                     %                     coeL_C(3,1) = WDpreA2*WignerD_single_element3*exp(1i*RightorLeft*mB*gamma)*coeL_B(1,1);
            %                     %                     coeL_C(4,1) = WDpreA2*WignerD_single_element4*exp(1i*RightorLeft*mB*gamma)*coeL_B(2,1);
            %                     %                     % for Y_l__m
            %                     %                     %coeL_C(1,1) = conj(coeL_C(1,1));
            %                     %                     %coeL_C(2,1) = conj(coeL_C(2,1));
            %                     %                     %coeL_C(3,1) = conj(coeL_C(3,1));
            %                     %                     %coeL_C(4,1) = conj(coeL_C(4,1));
            %                     %                     coeL_C = conj(coeL_C);
            %                     %                     if immproper
            %                     %                         %coeL_C(1,1) = (-1)^(A.l)*coeL_C(1,1);
            %                     %                         %coeL_C(2,1) = (-1)^(A.l)*coeL_C(2,1);
            %                     %                         coeL_C =  (-1)^(A.l)*coeL_C;
            %                     %                     end
            %                     %                     coeL = [coeL;coeL_C];
            %                     % back to realSH
            %                     if mB < 0
            %                         coeL_B(1,1)  =  1/sqrt(2);
            %                         coeL_B(2,1)  =  -1i/sqrt(2);
            %                     elseif mB == 0
            %                         coeL_B(1,1)  =  1/2;
            %                         coeL_B(2,1)  =  1/2;
            %                     elseif mB > 0
            %                         coeL_B(1,1)  =  1/sqrt(2)*(-1)^mB;
            %                         coeL_B(2,1)  =  1i/sqrt(2)*(-1)^mB;
            %                     end
            %                     mL = [mL;abs(mB);-abs(mB);abs(mB);-abs(mB);];
            %                     WignerD_single_element1 = (Y_l__m.d(Al,-A.m,mB,beta));
            %                     WignerD_single_element2 = (Y_l__m.d(Al, A.m,mB,beta));
            %                     coeL_C(1,1) = WDpreA1*WignerD_single_element1*exp(1i*RightorLeft*mB*gamma);
            %                     coeL_C(2,1) = WDpreA2*WignerD_single_element2*exp(1i*RightorLeft*mB*gamma);
            %                     % for Y_l__m
            %                     coeL_C(1,1) = conj(coeL_C(1,1));
            %                     coeL_C(2,1) = conj(coeL_C(2,1));
            %                     if immproper
            %                         coeL_C(1,1) = (-1)^(A.l)*coeL_C(1,1);
            %                         coeL_C(2,1) = (-1)^(A.l)*coeL_C(2,1);
            %                     end
            %                     coeL = [coeL;kron(coeL_C,coeL_B)];
            %                 end
            %                 [mL,coeL] = HollowKnight.generalcontractrow(mL,coeL);
            %             end

            %%%%%%%%%%% complex to real SH finish %%%%%%%%%%% 
            % kron
            nszL = size(szL,1);
            nmL = size(mL,1);
            mL = kron(ones(nszL,1),mL);
            coeL = kron(ones(nszL,1),coeL);
            szL = kron(szL,ones(nmL,1));
            scoeL = kron(scoeL,ones(nmL,1));
            %
            if ~conjugate
                CoeL = Coe_global*coeL.*scoeL;
            else
                CoeL = conj(Coe_global*coeL.*scoeL);
            end
            %% 
            QmumObj = Qnum.QnumL(A.AN,A.n,Al,mL,As,szL,A.q,CoeL,'nQnums',1,'nQnums_col',length(mL));
            %
            %disp(QmumObj);
        end
    end
    %%
    methods % disp
        function Qnumstring = string(Qnums)
            Qnumstring = mat2str([Qnums.AN,Qnums.l,Qnums.m,Qnums.s,Qnums.sz]);
        end
    end
    %% overload
    methods
        function C = eq(A,B,options)
            arguments
                A Qnum
                B Qnum
                options.strict = false
            end
            C = true;
            if A.AN ~= B.AN
                C =false;
                return;
            end
            if A.l ~= B.l
                C =false;
                return;
            end
            if A.m ~= B.m
                C =false;
                return;
            end
            if A.s ~= B.s
                C =false;
                return;
            end
            if A.sz ~= B.sz
                C =false;
                return;
            end
            if A.q ~= B.q
                C =false;
            end
            %
            if options.strict
                if A.coe ~= B.coe
                    C =false;
                end
            end
        end
    end
    %% contract
    methods
        function B = contractrow(A,options)
            arguments
                A Qnum;
                options.forgetcoe = false;
            end
            %
            B = A;%keep origin
            %cleanzero first
            A = cleanrow(A);
            % delete zero coe!!
            % ANL,nL,lL,mL,sL,szL,qL
            QnumL = ([A(:).AN;A(:).n;A(:).l;A(:).m;A(:).s;A(:).sz;A(:).q]).';
            coeL = ([A(:).coe]).';
            [QnumL_unique,coeL_unique] = HollowKnight.generalcontractrow(QnumL,coeL);
            [B,~,~] = Qnum.QnumL(...
                QnumL_unique(:,1),...
                QnumL_unique(:,2),...
                QnumL_unique(:,3),...
                QnumL_unique(:,4),...
                QnumL_unique(:,5),...
                QnumL_unique(:,6),...
                QnumL_unique(:,7),...
                coeL_unique,'nQnums',1,'nQnums_col',size(QnumL_unique,1));
            B = cleanrow(B);
        end
    end
    methods %get
        function hollow = get.hollow(QnumObj)
            if isnan(QnumObj.coe)
                hollow = true;
            else
                hollow = false;
            end
        end
    end
    methods %math
        
    end
end

