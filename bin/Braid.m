classdef Braid < matlab.mixin.CustomDisplay
    %UNTITLED7 Summary of this class goes here 
    %   Detailed explanation goes here

    properties
        BraidWord;
        Nstrings;
        Crossings;
        NLink;
    end
    properties
        Generator;
        BraidNum;
        ChP;
        Kesi_k;
        Hsym;
        Hvasplib;
    end
    properties
        NaiveSeperateBands;
        NaiveEIGENCAR;
        PlotSeqL;
        PlotSeqLMirror;
    end
    properties
        Permutation;
        CycleDecompositionStr;
        CycleDecomposition;
        CycleLift;
    end
    properties
        %CheckSignFunction;
    end
    properties % Data
        D_C = cell(1);
        kD_C = cell(1);
        ColorD_C =cell(1);
        Fnk = cell(1);
        Fjnk;
        kcross;
        Fkcross;
        kp;
        Gkp;
        Gnk = cell(1);
        Gjnk;
    end

    %% Define which properties show
    methods (Access = protected)
        function propgrp = getPropertyGroups(~)
            proplist = {'BraidWord','Nstrings','Crossings','NLink','CycleDecompositionStr'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    methods
        function BraidObj = Braid(BraidWord,Nstrings)
            DoubleBraidWord_col = Braid.BraidWord2Num(BraidWord);
            if nargin < 2
                Nstrings = Braid.NumBraidWord2Nstrings(DoubleBraidWord_col);
            end
            Crossings = length(DoubleBraidWord_col);
            [Generator,GeneratorStr] = Braid.Nstrings2Generator(Nstrings);
            Permutation =  Braid.PermutationBraidNum(DoubleBraidWord_col);
            [CycleDecomposition,CyclePresentationStr,CycleLift] = Braid.Cauchy2Cycle(Permutation);
            %
            BraidObj.BraidWord = BraidWord;
            BraidObj.BraidNum = DoubleBraidWord_col;
            BraidObj.Nstrings = Nstrings;
            BraidObj.Crossings = Crossings;
            BraidObj.Generator = Generator;
            BraidObj.Permutation = Permutation;
            BraidObj.CycleDecomposition = CycleDecomposition;
            BraidObj.CycleLift = CycleLift;
            BraidObj.NLink = length(CycleDecomposition);
            BraidObj.CycleDecompositionStr = CyclePresentationStr;
            [BraidObj.NaiveSeperateBands,BraidObj.NaiveEIGENCAR,BraidObj.PlotSeqL,BraidObj.PlotSeqLMirror] = NaiveBands(BraidObj);
        end

    end
    %%
    methods 
        function [Sign,k_index] = CheckCrossSign(BraidObj,kpoint)
            BraidNum = BraidObj.BraidNum;
            Cross = length(BraidNum); 
            BraidNumsign = sign(BraidNum);
            k_normal = kpoint./(2*pi) * Cross;
            k_index  = ceil(k_normal);
            Sign = BraidNumsign(k_index);
        end

    end
    %% Parent Function
    methods
        function BraidObj = D_C_gen(BraidObj)
            NormalPbands = BraidObj.NaiveEIGENCAR;
            Width = (max(NormalPbands,[],"all")-min(NormalPbands,[],"all"))/2;
            Middle = (max(NormalPbands,[],"all")+ min(NormalPbands,[],"all"))/2;
            NormalPbands = (-1)*(NormalPbands-Middle)/Width;
            Lk = BraidObj.CycleDecomposition;
            CycleLift = BraidObj.CycleLift;
            Nsn = length(Lk);
            %CK = BraidObj.Crossings;
            D_Cj = NormalPbands;
            %l = CK;
            %Listk = 1:l;
            %t_k = 2*sym(pi)*(2*Listk - 1)/2*l;
            %PermutationTmp = BraidObj.Permutation(2,:);
            % as DOI: 10.1142/S0218216518500827 Convention:
            % 3. Finding Fourier Parametrizations
            % Finding the data points for the trigonometric interpolation for FC
            for n = 1:Nsn
                CycleTmpLift = CycleLift{n};
                s = length(CycleTmpLift);
                D_C{n} = [D_Cj(CycleTmpLift(1),:)];
                ColorD_C{n} = [CycleTmpLift(1)];
                %
                for j = 2:s
                    D_C{n} = [D_C{n},D_Cj(CycleTmpLift(j),2:end);];
                    ColorD_C{n} = [ColorD_C{n},CycleTmpLift(j)];
                end
                kD_C{n} = linspace(0,2*sym(pi),size(D_C{n},2));
                %kD_C{n}(end) = false;
            end
            BraidObj.D_C = D_C;
            BraidObj.kD_C = kD_C;
            BraidObj.ColorD_C = ColorD_C;
        end
        function [BraidObj,am_n] = FnkGen(BraidObj)
            if isempty(BraidObj.Fjnk)
                BraidObj = D_C_gen(BraidObj);
            end
            Lk = BraidObj.CycleDecomposition;
            Cross = BraidObj.Crossings;
            for n = 1:numel(BraidObj.D_C)
                ln = length(Lk{n});
                Crossln = Cross*ln;
                kD_Cn = BraidObj.kD_C{n}(1:end-1);
                D_Cn = BraidObj.D_C{n}(1:end-1);
                if  mod(Crossln,2) == 1
                    mL = (- Crossln/2 + 1/2 ) : ( Crossln/2 - 1/2 ) ;
                else
                    mL = (- Crossln/2 + 1) : ( Crossln/2 - 1 ) ;
                    m = Crossln/2;
                    expL =  exp(-1i.*kD_Cn*m);
                    a_Crossln_2 = 1/(Crossln) *  sum(D_Cn.*expL);
                end
                am_n{n} = zeros(1,length(mL));
                count = 0;
                for m = mL
                    count = count + 1;
                    expL =  exp(-1i.*kD_Cn*m);
                    am = 1/(Crossln) *  sum(D_Cn.*expL);
                    am_n{n}(count) = am;
                end
                syms k real;
                expkL = exp(1i*mL*k);
                if mod(Crossln,2) == 1
                    Fnk{n} = simplify(rewrite(sum(expkL.*am_n{n}),'sincos'));
                else
                    Fnk{n} = simplify(rewrite(sum([expkL.*am_n{n},a_Crossln_2*cos(Crossln/2 * k)]),"sincos"));
                end
                kjn = (k+2*sym(pi)*(0:(ln-1)))/ln;
                for ikj = 1:numel(kjn)
                    Fjnk{n}{ikj} = subs(Fnk{n},k,kjn(ikj));
                end
            end
            BraidObj.Fnk = Fnk;
            BraidObj.Fjnk = Fjnk;
        end
        function BraidObj = kpGen(BraidObj,kn)
            arguments
                BraidObj
                kn = 100;
            end
            if isempty(BraidObj.kcross)
                BraidObj = FnkGen(BraidObj);
            end
            syms k real;
            kL  = linspace(0,2*pi,kn);
            Lk = BraidObj.CycleDecomposition;
            CycleLift = BraidObj.CycleLift;
            NaiveEIGENCAR = BraidObj.NaiveEIGENCAR;
            %Cross = BraidObj.Crossings;
            %kCross = linspace(0,2*sym(pi),Cross + 2);
            %kCrossInterval = 
            %kCross([1,end]) = [];BraidObj.nStrings
            %
            Falljn_label = [];
            count = 0;
            for n = 1:numel(BraidObj.Fjnk)
                ln = length(BraidObj.Fjnk{n});
                for j = 1:ln
                    count = count + 1;
                    Falljn{count} = BraidObj.Fjnk{n}{j};
                    Falljn_label = [Falljn_label;[n,j,ln]];
                end
            end
            %
            kcross{n} = [];
            Fkcross{n} = [];
            kp{n} = [];
            Gkp{n} = [];
            nF = count;
            %
            for i = 1:nF
                TargetFjnk_left = Falljn{i};
                n = Falljn_label(i,1);
                CycleLiftn = CycleLift{n};
                jforstring = Falljn_label(i,2);
                iString = CycleLiftn(jforstring);
                n_using = n;
                ln_using = Falljn_label(i,3);
                j_using = jforstring;
                for j = 1:nF
                    if i == j
                        continue;
                    end
                    TargetFjnk_right = Falljn{j};
                    n = Falljn_label(j,1);
                    CycleLiftn = CycleLift{n};
                    jforstring = Falljn_label(j,2);
                    jString = CycleLiftn(jforstring);
                    %
                    EQtmpZero = simplify(TargetFjnk_left - TargetFjnk_right);
                    EQtmp = EQtmpZero == 0;
                    % numerical
                    EQtmpF = matlabFunction(EQtmpZero);
                    EQtmpL = EQtmpF(kL);
                    EQtmpsignL = sign(EQtmpL);
                    EQtmpsigndiffL = diff(EQtmpsignL);
                    CrossLocation = find(EQtmpsigndiffL);
                    if ~isempty(CrossLocation)
                        for d = CrossLocation
                            syms k real;
                            EQtmp2 = (k>sym(kL(d-1)));
                            EQtmp3 = (k<sym(kL(d+1)));
                            kcrosstmp = solve([EQtmp,EQtmp2,EQtmp3],k);
                            try
                                double(kcrosstmp);
                            catch 
                                kcrosstmp = [];
                            end
                            if isempty(kcrosstmp)
                                kcrosstmp = vpasolve(EQtmp,k,[kL(d-1),kL(d+1)]);
                            end
                            Fkcrosstmp = subs(TargetFjnk_left,k,kcrosstmp);
                            [TmpSign,k_index] =  BraidObj.CheckCrossSign(kcrosstmp);
                            StringSeq = NaiveEIGENCAR([iString,jString],k_index);
                            StringSeqSign = sign(StringSeq(2)-StringSeq(1));
                            kp_left = (kcrosstmp+2*pi*(j_using-1))/ln_using;
                            %kp_j = (kcrosstmp+2*pi*(j-1))/ln;
                            if TmpSign*StringSeqSign < 0
                                Gkp_left = 1;
                                %Gkp_j = -1;
                            else
                                Gkp_left = -1;
                                %Gkp_j = 1;
                            end
                            kcross{n_using} = [kcross{n_using},kcrosstmp];
                            Fkcross{n_using} = [Fkcross{n_using},Fkcrosstmp];
                            kp{n_using} = [kp{n_using},[kp_left ]];
                            Gkp{n_using} = [Gkp{n_using},[Gkp_left ]];
                        end
                    end
                end
            end
            for n = 1:numel(BraidObj.Fjnk)
                [~,order] =  sort(double(kcross{n}));
                Fkcross{n} = Fkcross{n}(order);
                kcross{n} = kcross{n}(order);
                [~,order2] =  sort(double(kp{n}));
                kp{n} = kp{n}(order2);
                Gkp{n} = Gkp{n}(order2);
            end
            % for n = 1:numel(BraidObj.Fjnk)
            %     CycleLiftn = CycleLift{n};
            %     TargetFjnk = BraidObj.Fjnk{n};
            %     ln = length(TargetFjnk);
            %     kcross{n} = [];
            %     Fkcross{n} = [];
            %     kp{n} = [];
            %     Gkp{n} = [];
            %     % try symbolick
            %     for i = 1:ln
            %         iString = CycleLiftn(i);
            %         for j = i:ln
            %             jString = CycleLiftn(j);
            %             if iString == jString
            %                 continue;
            %             end
            %             EQtmpZero = simplify(TargetFjnk{i} - TargetFjnk{j});
            %             EQtmp = EQtmpZero == 0;
            %             % numerical
            %             EQtmpF = matlabFunction(EQtmpZero);
            %             EQtmpL = EQtmpF(kL);
            %             EQtmpsignL = sign(EQtmpL);
            %             EQtmpsigndiffL = diff(EQtmpsignL);
            %             CrossLocation = find(EQtmpsigndiffL);
            %             if ~isempty(CrossLocation)
            %                 for d = CrossLocation
            %                     syms k real;
            %                     EQtmp2 = (k>sym(kL(d-1)));
            %                     EQtmp3 = (k<sym(kL(d+1)));
            %                     kcrosstmp = solve([EQtmp,EQtmp2,EQtmp3],k);
            %                     if isempty(kcrosstmp)
            %                         kcrosstmp = vpasolve([EQtmp],k,[kL(d-1),kL(d+1)]);
            %                     end
            %                     Fkcrosstmp = subs(TargetFjnk{i},k,kcrosstmp);
            %                     [TmpSign,k_index] =  BraidObj.CheckCrossSign(kcrosstmp);
            %                     StringSeq = NaiveEIGENCAR([iString,jString],k_index);
            %                     StringSeqSign = sign(StringSeq(2)-StringSeq(1));
            %                     kp_i = (kcrosstmp+2*pi*(i-1))/ln;
            %                     kp_j = (kcrosstmp+2*pi*(j-1))/ln;
            %                     if TmpSign*StringSeqSign < 0
            %                         Gkp_i = 1;
            %                         Gkp_j = -1;
            %                     else
            %                         Gkp_i = -1;
            %                         Gkp_j = 1;
            %                     end
            %                     kcross{n} = [kcross{n},kcrosstmp];
            %                     Fkcross{n} = [Fkcross{n},Fkcrosstmp];
            %                     kp{n} = [kp{n},[kp_i kp_j]];
            %                     Gkp{n} = [Gkp{n},[Gkp_i Gkp_j]];
            %                 end
            %             end
            %         end
            %     end
            %     [~,order] =  sort(double(kcross{n}));
            %     Fkcross{n} = Fkcross{n}(order);
            %     kcross{n} = kcross{n}(order);
            %     [~,order2] =  sort(double(kp{n}));
            %     kp{n} = kp{n}(order2);
            %     Gkp{n} = Gkp{n}(order2);
            % end
            BraidObj.kcross = kcross;
            BraidObj.Fkcross = Fkcross;
            BraidObj.kp = kp;
            BraidObj.Gkp = Gkp;
        end
        function BraidObj = GkData_gen(BraidObj)
            kpset = BraidObj.kcross;
            kpsetforGk = kpset;
            syms k real;
            for n = 1:numel(kpset)
                ln = length(BraidObj.CycleDecomposition{n});
                tmpkset = [];
                for j = 0:ln-1
                    tmpkset = [tmpkset,(kpset{n}+2*pi*j)/ln];
                end
                for i = 1:numel(tmpkset)
                    Gkpset{n}(i) = BraidObj.CheckCrossSign(tmpkset(i));
                end
                kpsetforGk{n} = tmpkset;
            end
            BraidObj.kpsetforGk = kpsetforGk;
            BraidObj.Gkpset = Gkpset;
        end
        function [BraidObj,bm_n] = GnkGen2(BraidObj)
            % The Method from PHYSICAL REVIEW LETTERS 126, 010401 (2021)
            % uses X = A\b;
            % Try 10.1142/S0218216518500827 Convention:
            % Step 4: Trigonometric interpolation for GC
            % The reason is also discussied
            % "Since the positions of the crossings of B′ are in general not
            % equidistributed, the trigonometric interpolation does not
            % directly translate to a discrete Fourier Transform."
            % we use j -> k
            % we use k -> t
            if isempty(BraidObj.Gkp)
                BraidObj = kpGen(BraidObj);
            end
            Lk = BraidObj.CycleDecomposition;
            for n = 1:numel(BraidObj.Gkp)
                yk = BraidObj.Gkp{n};
                ln = length(Lk{n});
                N = length(yk);
                tkL = BraidObj.kp{n};
                Odd = mod(N,2) == 1;
                jL =  1:N;
                syms k real;
                expL = exp(1i.*tkL);
                ChooseL = 1:N;
                %Cumprod = 
                if Odd
                    K = (N-1)/2;
                    CoeffL = yk(jL) .* exp((-1i * K * k) + (1i * K * tkL));
                    count = 0;
                    for tk_prime = tkL
                        count = count + 1;
                        ChooseLTmp  =ChooseL;
                        ChooseLTmp(ChooseL == count) = [];
                        expLmodify = expL(ChooseLTmp);
                        Cumprod(count) = fold(@times,(exp(1i*k) - expLmodify)./(exp(1i*tk_prime) - expLmodify));
                    end
                    Gnk{n} = simplify(sum(CoeffL.*Cumprod));
                else
                    K = (N)/2;
                    CoeffL = yk(jL) .* exp((-1i * K * k) + (1i * K * tkL));
                    count = 0;

                    for tk_prime = tkL
                        count = count + 1;
                        ChooseLTmp  =ChooseL;
                        ChooseLTmp(ChooseL == count) = [];
                        expLmodify = expL(ChooseLTmp);
                        expLmodify = [1,expLmodify];
                        Cumprod(count) = fold(@times,(exp(1i*k) - expLmodify)./(exp(1i*tk_prime) - expLmodify));
                    end
                    Gnk{n} = simplify(sum(CoeffL.*Cumprod));
                end
                kjn = (k+2*sym(pi)*(0:(ln-1)))/ln;
                for ikj = 1:numel(kjn)
                    Gjnk{n}{ikj} = subs(Gnk{n},k,kjn(ikj));
                end
            end
            BraidObj.Gnk = Gnk;
            BraidObj.Gjnk = Gjnk;
        end
        function [BraidObj,bm_n] = GnkGen(BraidObj)
            % The Method from PHYSICAL REVIEW LETTERS 126, 010401 (2021)
            if isempty(BraidObj.Gkp)
                BraidObj = kpGen(BraidObj);
            end
            Lk = BraidObj.CycleDecomposition;
            for n = 1:numel(BraidObj.Gkp)
                D_Cn = BraidObj.Gkp{n};
                ln = length(Lk{n});
                N = length(D_Cn);
                kD_Cn = BraidObj.kp{n};
                Odd = mod(N,2) == 1;
                %jL =  1:N;
                syms k real;
                %expL = exp(1i.*kD_Cn);
                Crossln = N;
                Crossln_2 = N/2;
                if  Odd
                    mL = (- Crossln/2 + 1/2 ) : ( Crossln/2 - 1/2 ) ;
                    Width = length(mL);
                    %bmTosolve = sym('b_m',[1 Width]);
                    expSymL =   [exp(1i*k.*mL)];
                    expLFunction = matlabFunction(expSymL,'Vars',k);
                else
                    mL = (- Crossln/2 + 1) : ( Crossln/2 - 1 ) ;
                    Width = length(mL) + 1;
                    %bmTosolve = sym('b_m',[1 Width]);
                    expSymL =   [exp(1i*k.*mL) cos(Crossln_2*k)];
                    expLFunction = matlabFunction(expSymL,'Vars',k);
                end
                TheMat = zeros([N,Width],'sym');
                for i = 1:N
                    TheMat(i,:) =  expLFunction(kD_Cn(i));
                end
                % Try Cramer
                TheMat = simplify(TheMat);
                if double(det(TheMat)) == 0
                    yk = BraidObj.Gkp{n};
                    ln = length(Lk{n});
                    tkL = BraidObj.kp{n};
                    jL =  1:N;
                    expL = exp(1i.*tkL);
                    ChooseL = 1:N;
                    if Odd
                        K = (N-1)/2;
                        CoeffL = yk(jL) .* exp((-1i * K * k) + (1i * K * tkL));
                        count = 0;
                        for tk_prime = tkL
                            count = count + 1;
                            ChooseLTmp  =ChooseL;
                            ChooseLTmp(ChooseL == count) = [];
                            expLmodify = expL(ChooseLTmp);
                            Cumprod(count) = fold(@times,(exp(1i*k) - expLmodify)./(exp(1i*tk_prime) - expLmodify));
                        end
                        Gnk{n} = simplify(sum(CoeffL.*Cumprod));
                    else
                        K = (N)/2;
                        CoeffL = yk(jL) .* exp((-1i * K * k) + (1i * K * tkL));
                        count = 0;

                        for tk_prime = tkL
                            count = count + 1;
                            ChooseLTmp  =ChooseL;
                            ChooseLTmp(ChooseL == count) = [];
                            expLmodify = expL(ChooseLTmp);
                            expLmodify = [1,expLmodify];
                            Cumprod(count) = fold(@times,(exp(1i*k) - expLmodify)./(exp(1i*tk_prime) - expLmodify));
                        end
                        Gnk{n} = simplify(sum(CoeffL.*Cumprod));
                    end
                else
                    TheSolve = (TheMat)\D_Cn.';
                    Gnk{n} = simplify(sum((TheSolve.').*expSymL));
                end

                kjn = (k+2*sym(pi)*(0:(ln-1)))/ln;
                for ikj = 1:numel(kjn)
                    Gjnk{n}{ikj} = subs(Gnk{n},k,kjn(ikj));
                end
            end
            BraidObj.Gnk = Gnk;
            BraidObj.Gjnk = Gjnk;
        end

        function BraidObj = ChPGen(BraidObj)
            syms lambda;
            syms k real;
            Flambdak = sym(1);
            if isempty(BraidObj.Gjnk)
                BraidObj = BraidObj.GnkGen();
            end
            %NLink = length(BraidObj.NLink); 
            FjnCell = BraidObj.Fjnk;
            GjnCell = BraidObj.Gjnk;
            for n = 1:BraidObj.NLink
                ln = length(BraidObj.CycleDecomposition{n});
                for j = 1:ln
                    Flambdak = Flambdak*(lambda-FjnCell{n}{j}-1i*GjnCell{n}{j});
                end
            end
            Flambdak = simplify(Flambdak);
            BraidObj.ChP = Flambdak;
            LaurentSeries = series(BraidObj.ChP,lambda,'Order',BraidObj.Nstrings);
            BraidObj.Kesi_k  = coeffs(LaurentSeries,lambda);
            %= flip(C);
        end
        function BraidObj = HsymGen(BraidObj)
            syms k k_x real;
            if isempty(BraidObj.ChP)
                BraidObj = BraidObj.ChPGen();
            end
            Hsym = sym(diag(ones(1,BraidObj.Nstrings-1),-1));
            Hsym(1,:) = -flip( BraidObj.Kesi_k);
            BraidObj.Hsym = subs(Hsym,k,k_x);
        end
        function BraidObj = vasplibGen(BraidObj)
            syms k k_x real;
            Hvasplib = vasplib();
            Hvasplib.Basis_num = BraidObj.Nstrings;
            Hvasplib.Rm = [1];
            Hvasplib.Dim = 1;
            Hvasplib.Hsym = BraidObj.Hsym;
            kpoints_cart = [
                0;
                pi;

                pi;
                2*pi;
                ];
            kpoints_frac = kpoints_cart/(2*pi);

            [Hvasplib.klist_cart,Hvasplib.klist_frac,Hvasplib.klist_l,Hvasplib.kpoints_l,~] = ...
                vasplib.kpathgen(kpoints_frac, ...
                61,Hvasplib.Gk,'Dim',1);
            Hvasplib.kpoints_name = ["0","\pi","2\pi"];
            BraidObj.Hvasplib = Hvasplib;
        end
    end
    %% DataSet 
    methods
        function [NaiveSeperateBands,NaiveEIGENCAR,PlotSeqL,PlotSeqLVertical] = NaiveBands(BraidObj)
            NaiveSeperateBands{BraidObj.Crossings} = (1:BraidObj.Nstrings).';
            PermutationList = abs(BraidObj.BraidNum);
            BraidList = BraidObj.BraidNum;
            Nbands = BraidObj.Nstrings;
            Ncross = BraidObj.Crossings;
            DefaultL = (1:BraidObj.Nstrings).';
            PlotSeqL = repmat(DefaultL,[1 Ncross]);
            PlotSeqLVertical = PlotSeqL;
            NaiveEIGENCAR = zeros(Nbands,Ncross+1);
            NaiveEIGENCAR(:,1) = DefaultL;
            for i = 1:Ncross
                if i > 1
                    NaiveSeperateBands{i}(:,1) =  NaiveSeperateBands{i-1}(:,2);
                else
                    NaiveSeperateBands{i}(:,1) = DefaultL;
                end
                PlusN = -1;
                MinusN = -1;
                for n = 1:Nbands
                    if NaiveSeperateBands{i}(n,1) == PermutationList(i)
                        NaiveSeperateBands{i}(n,2) = NaiveSeperateBands{i}(n,1) + 1;
                        PlusN = n;
                    elseif NaiveSeperateBands{i}(n,1) == PermutationList(i) +1
                        NaiveSeperateBands{i}(n,2) = NaiveSeperateBands{i}(n,1) - 1;
                        MinusN = n;
                    else
                        NaiveSeperateBands{i}(n,2) = NaiveSeperateBands{i}(n,1);
                    end
                end
                if PlusN > 0 && MinusN > 0
                    if BraidList(i) > 0 && PlusN > MinusN  
                        
                    elseif  BraidList(i) < 0 && PlusN < MinusN 
                    else
                        PlotSeqL([MinusN,PlusN],i) = ...
                            PlotSeqL([PlusN,MinusN],i);
                        PlotSeqLVertical([MinusN,PlusN],i)= ...
                            PlotSeqLVertical([PlusN,MinusN],i);
                    end
                end
                %BeginBand = NaiveSeperateBands{i}(:,1);
                %PlotSeqL(:,i) =  NaiveSeperateBands{i}(:,1);
                %PlotSeqLVertical(:,i)  = NaiveSeperateBands{i}(:,1);
                NaiveEIGENCAR(:,i+1) = NaiveSeperateBands{i}(:,2);
            end
        end
    end
    %% plot
    methods
        function ax = Fjnk_Plot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                %options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 10;
                options.MarkerSize = -1;
                options.MarkerEdgeColor = [1/2,1/2,1/2];
                options.MarkerFaceColor = [1/2,1/2,1/2];
                opts.vertical = false;
                opts.kn = 100;
            end
            Nbands=BraidObj.Nstrings;
            %BraidObj = D_C_gen(BraidObj);
            if isempty(BraidObj.Fkcross)
                BraidObj = BraidObj.kpGen();
            end
            D_C = BraidObj.D_C;
            kD_C = BraidObj.kD_C;
            ColorD_C = BraidObj.ColorD_C;
            Fjnk = BraidObj.Fjnk;
            [~,Ax] = Figs(length(BraidObj.D_C),1);
            %ax = options.ax;
            %hold(ax,'on');
            if options.MarkerSize <0
                options.MarkerSize  = options.LineWidth*50;
            end
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for n = 1:Nbands
                    colorL(n,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);

            if opts.vertical
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = size(D_C{n},2)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = latex(kD_C{n}(i));
                    end
                    for j = 1:length(ColorD_C{n})
                        FjnkFunction = matlabFunction(Fjnk{n}{j},'Vars',sym('k'));
                        Ei = ColorD_C{n}(j);
                        KPs = [0,2*pi];
                        KL = linspace(KPs(1),KPs(end),opts.kn);
                        FL = FjnkFunction(KL);
                        plot(ax,double(FL),double(KL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'DisplayName',['F_',num2str(n),'^',num2str(j),'(k)']);
                        SelectL = SelectL + Ndivide -1;
                    end
                    scatter(ax,...
                        double(BraidObj.Fkcross{n}),...
                        double(BraidObj.kcross{n}),...
                        options.MarkerSize,...
                        options.MarkerFaceColor,...
                        'filled','DisplayName','Cross');
                    ylabel(ax,'k');
                    yticks(ax,ticks);
                    set(ax,'TickLabelInterpreter','latex');
                    set(ax,'FontName','Times');
                    yticklabels(ax,ticklabel);
                    ylim(ax,[0-pi/10,2*pi+pi/10]);
                    xlim(ax,[-1.2,1.2]);
                    xlabel(ax,'x');
                    xticks(ax,[-1,0,1]);
                end
            else
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = 1+ (size(D_C{n},2)-1)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    set(ax,'FontName','Times');
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = "$"+latex(kD_C{n}(i))+"$";
                    end
                    for j = 1:length(ColorD_C{n})
                        FjnkFunction = matlabFunction(Fjnk{n}{j},'Vars',sym('k'));
                        Ei = ColorD_C{n}(j);
                        KPs = [0,2*pi];
                        KL = linspace(KPs(1),KPs(end),opts.kn);
                        FL = FjnkFunction(KL);
                        plot(ax,double(KL),double(FL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'DisplayName',['F_',num2str(n),'^',num2str(j),'(k)']);
                        %fplot(ax,[0,2*pi],options.LineSpec,"LineWidth",options.LineWidth);
                        SelectL = SelectL + Ndivide -1;
                    end
                    scatter(ax,double(BraidObj.kcross{n}),...
                        double(BraidObj.Fkcross{n}),...
                        options.MarkerSize,...
                        options.MarkerFaceColor,...
                        'filled','DisplayName','Cross');
                    set(ax,'TickLabelInterpreter','latex');
                    ylabel(ax,'y');
                    yticks(ax,[-1,0,1]);
                    xlabel(ax,'k');
                    xticks(ax,ticks);
                    xlim(ax,[0-pi/10,2*pi+pi/10]);
                    ylim(ax,[-1.2,1.2]);
                    xticklabels(ax,ticklabel);
                end
            end

            %axis(ax,'off');
        end
        function ax = Fnk_Plot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                %options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 10;
                options.MarkerSize = -1;
                options.MarkerEdgeColor = [1/2,1/2,1/2];
                options.MarkerFaceColor = 'none';
                opts.vertical = false;
                opts.kn = 100;
            end
            Nbands=BraidObj.Nstrings;
            %BraidObj = D_C_gen(BraidObj);
            D_C = BraidObj.D_C;
            kD_C = BraidObj.kD_C;
            ColorD_C = BraidObj.ColorD_C;
            Fnk = BraidObj.Fnk;
            [~,Ax] = Figs(length(BraidObj.D_C),1);
            %ax = options.ax;
            %hold(ax,'on');
            if options.MarkerSize <0
                options.MarkerSize  = options.LineWidth*3;
            end
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for n = 1:Nbands
                    colorL(n,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);

            if opts.vertical
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = size(D_C{n},2)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = latex(kD_C{n}(i));
                    end
                    FnkFunction = matlabFunction(Fnk{n},'Vars',sym('k'));
                    for j = 1:length(ColorD_C{n})
                        Ei = ColorD_C{n}(j);
                        KPs = kD_C{n}(SelectL);
                        KL = linspace(KPs(1),KPs(end),opts.kn);
                        FL = FnkFunction(KL);
                        plot(ax,double(FL),double(KL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'DisplayName',['EFit_',num2str(Ei)]);
                        plot(ax,double(D_C{n}(SelectL)),double(KPs),'o',...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',ModifycolorL(Ei,:),...
                            'DisplayName',['E_',num2str(Ei)]);
                        SelectL = SelectL + Ndivide -1;
                    end
                    ylabel(ax,'k');
                    yticks(ax,ticks);
                    set(ax,'TickLabelInterpreter','latex');
                    set(ax,'FontName','Times');
                    yticklabels(ax,ticklabel);
                    ylim(ax,[0-pi/10,2*pi+pi/10]);
                    xlim(ax,[-1.2,1.2]);
                    xlabel(ax,'x');
                    xticks(ax,[-1,0,1]);
                end
            else
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = 1+ (size(D_C{n},2)-1)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    set(ax,'FontName','Times');
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = "$"+latex(kD_C{n}(i))+"$";
                    end
                    FnkFunction = matlabFunction(Fnk{n},'Vars',sym('k'));
                    for j = 1:length(ColorD_C{n})
                        Ei = ColorD_C{n}(j);
                        KPs = kD_C{n}(SelectL);
                        KL = linspace(KPs(1),KPs(end),opts.kn);
                        FL = FnkFunction(KL);
                        plot(ax,double(KL),double(FL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'DisplayName',['EFit_',num2str(Ei)]);
                        plot(ax,double(KPs),D_C{n}(SelectL),'o',...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',ModifycolorL(Ei,:),...
                            'DisplayName',['E_',num2str(Ei)]);
                        %fplot(ax,[0,2*pi],options.LineSpec,"LineWidth",options.LineWidth);
                        SelectL = SelectL + Ndivide -1;
                    end
                    set(ax,'TickLabelInterpreter','latex');
                    ylabel(ax,'y');
                    yticks(ax,[-1,0,1]);
                    xlabel(ax,'k');
                    xticks(ax,ticks);
                    xlim(ax,[0-pi/10,2*pi+pi/10]);
                    ylim(ax,[-1.2,1.2]);
                    xticklabels(ax,ticklabel);
                end
            end

            %axis(ax,'off');
        end
        function ax = Gnk_Plot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                %options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 5;
                options.MarkerSize = -1;
                options.MarkerEdgeColor = [1/2,1/2,1/2];
                options.MarkerFaceColor = 'r';
                opts.vertical = false;
                opts.kn = 361;
            end
            Nbands=BraidObj.Nstrings;
            %BraidObj = D_C_gen(BraidObj);
            Gkp = BraidObj.Gkp;
            kp = BraidObj.kp;
            Gnk = BraidObj.Gnk;
            [~,Ax] = Figs(length(BraidObj.Gkp),1);
            %ax = options.ax;
            %hold(ax,'on');
            if options.MarkerSize <0
                options.MarkerSize  = options.LineWidth*3;
            end
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for n = 1:Nbands
                    colorL(n,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);

            if opts.vertical
                for n = 1:length(Ax)
                    ax = Ax(n);
                    ticks = double(kp{n});
                    ticklabel = repmat("",[1 length(kp{n})]);
                    for i = 1:length(kp{n})
                        ticklabel(i) = latex(kp{n}(i));
                    end
                    GnkFunction = matlabFunction(Gnk{n},'Vars',sym('k'));
                    KPs = [0,2*pi];
                    KL = linspace(KPs(1),KPs(end),opts.kn);
                    GL = GnkFunction(KL);
                    plot(ax,double(GL),double(KL),options.LineSpec,...
                        'LineWidth',options.LineWidth,...
                        'Color',options.MarkerEdgeColor,...
                        'DisplayName',['G_',num2str(n),'(k)']);
                    plot(ax,double(Gkp{n}),double(kp{n}),'o',...
                        'LineWidth',options.LineWidth,...
                        'Color',options.MarkerFaceColor ,...
                        'MarkerSize',options.MarkerSize,...
                        'MarkerEdgeColor',options.MarkerEdgeColor,...
                        'MarkerFaceColor',options.MarkerFaceColor ,...
                        'DisplayName',['G_',num2str(n),'(k_p)']);
                    ylabel(ax,'k');
                    yticks(ax,ticks);
                    set(ax,'TickLabelInterpreter','latex');
                    set(ax,'FontName','Times');
                    yticklabels(ax,ticklabel);
                    ylim(ax,[0-pi/10,2*pi+pi/10]);
                    %xlim(ax,[-1.2,1.2]);
                    xlabel(ax,'x');
                    xticks(ax,[-1,0,1]);
                end
            else
                for n = 1:length(Ax)
                    ax = Ax(n);
                    ticks = double(kp{n});
                    ticklabel = repmat("",[1 length(kp{n})]);
                    set(ax,'FontName','Times');
                    for i = 1:length(kp{n})
                        ticklabel(i) = "$"+latex(kp{n}(i))+"$";
                    end
                    GnkFunction = matlabFunction(Gnk{n},'Vars',sym('k'));
                    KPs = [0,2*pi];
                    KL = linspace(KPs(1),KPs(end),opts.kn);
                    GL = GnkFunction(KL);
                    plot(ax,double(KL),double(GL),options.LineSpec,...
                        'LineWidth',options.LineWidth,...
                        'Color',options.MarkerEdgeColor,...
                        'DisplayName',['G_',num2str(n),'(k)']);
                    plot(ax,double(kp{n}),double(Gkp{n}),'o',...
                        'LineWidth',options.LineWidth,...
                        'Color',options.MarkerFaceColor ,...
                        'MarkerSize',options.MarkerSize,...
                        'MarkerEdgeColor',options.MarkerEdgeColor,...
                        'MarkerFaceColor',options.MarkerFaceColor ,...
                        'DisplayName',['G_',num2str(n),'(k_p)']);
                    set(ax,'TickLabelInterpreter','latex');
                    ylabel(ax,'y');
                    yticks(ax,[-1,0,1]);
                    xlabel(ax,'k');
                    xticks(ax,ticks);
                    xlim(ax,[0-pi/10,2*pi+pi/10]);
                    %ylim(ax,[-1.2,1.2]);
                    xticklabels(ax,ticklabel);
                end
            end

            %axis(ax,'off');
        end
        
        function ax = D_C_Plot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                %options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-o';
                options.LineWidth = 10;
                options.MarkerSize = -1;
                options.MarkerEdgeColor = [1/2,1/2,1/2];
                options.MarkerFaceColor = 'none';
                opts.vertical = false;
            end
            Nbands=BraidObj.Nstrings;
            %BraidObj = D_C_gen(BraidObj);
            D_C = BraidObj.D_C;
            kD_C = BraidObj.kD_C;
            ColorD_C = BraidObj.ColorD_C;
            [~,Ax] = Figs(length(BraidObj.D_C),1);
            CycleLift = BraidObj.CycleLift;
            %ax = options.ax;
            %hold(ax,'on');
            if options.MarkerSize <0
                options.MarkerSize  = options.LineWidth*3;
            end
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for n = 1:Nbands
                    colorL(n,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);

            if opts.vertical
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = size(D_C{n},2)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    %CycleLiftn = CycleLift{n};
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = latex(kD_C{n}(i));
                    end
                    for j = 1:length(ColorD_C{n})
                        Ei = ColorD_C{n}(j);
                        plot(ax,D_C{n}(SelectL),kD_C{n}(SelectL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',ModifycolorL(Ei,:),...
                            'DisplayName',['E_',num2str(Ei)]);
                        SelectL = SelectL + Ndivide -1;
                    end
                    ylabel(ax,'k');
                    yticks(ax,ticks);
                    set(ax,'TickLabelInterpreter','latex');
                    set(ax,'FontName','Times');
                    yticklabels(ax,ticklabel);
                    ylim(ax,[0-pi/10,2*pi+pi/10]);
                    xlim(ax,[-1.2,1.2]);
                    xlabel(ax,'x');
                    xticks(ax,[-1,0,1]);
                end
            else
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = 1+ (size(D_C{n},2)-1)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    set(ax,'FontName','Times');
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = "$"+latex(kD_C{n}(i))+"$";
                    end
                    for j = 1:length(ColorD_C{n})
                        Ei = ColorD_C{n}(j);
                        plot(ax,double(kD_C{n}(SelectL)),D_C{n}(SelectL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',ModifycolorL(Ei,:),...
                            'DisplayName',['E_',num2str(Ei)]);
                        SelectL = SelectL + Ndivide -1;
                    end
                    set(ax,'TickLabelInterpreter','latex');
                    ylabel(ax,'y');
                    yticks(ax,[-1,0,1]);
                    xlabel(ax,'k');
                    xticks(ax,ticks);
                    xlim(ax,[0-pi/10,2*pi+pi/10]);
                    ylim(ax,[-1.2,1.2]);
                    xticklabels(ax,ticklabel);
                end
            end

            %axis(ax,'off');
        end
        function ax = LinePlot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 10;
                options.MarkerSize = 3;
                options.MarkerEdgeColor = 'none';
                options.MarkerFaceColor = 'none';
                opts.vertical = false;
            end
            Nbands=BraidObj.Nstrings;
            ax = options.ax;
            hold(ax,'on');
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for i = 1:Nbands
                    colorL(i,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);
            PlotSeqList = BraidObj.PlotSeqL;
            EIGENCARCELL = BraidObj.NaiveSeperateBands;
            EIGENCAR = BraidObj.NaiveEIGENCAR;
            Ncross = BraidObj.Crossings;       
            if opts.vertical
                PlotSeqList = BraidObj.PlotSeqL;
                for i = 1:Ncross
                    for n = 1:BraidObj.Nstrings
                        
                        Ei = PlotSeqList(n,i);
                        plot(ax,EIGENCARCELL{i}(Ei,:),[i-1 i],options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',options.MarkerFaceColor,...
                            'DisplayName',['E_',num2str(Ei)]);
                    end
                end
                for n = 1:Nbands
                    scatter(ax,EIGENCAR(n,:),0:Ncross,...
                        ones(1,Ncross+1)*options.LineWidth*20,ModifycolorL(n,:),"filled",'DisplayName',['String_',num2str(n)]);
                end
            else
                for i = 1:Ncross
                    EIGENCAR(:,i+1) = EIGENCARCELL{i}(:,2);
                    for n = 1:BraidObj.Nstrings
                        
                        Ei = PlotSeqList(n,i);
                        plot(ax,[i-1 i],-EIGENCARCELL{i}(Ei,:),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',options.MarkerFaceColor,...
                            'DisplayName',['E_',num2str(Ei)]);
                    end
                end
                for n = 1:Nbands
                    scatter(ax,0:Ncross,-EIGENCAR(n,:),...
                        ones(1,Ncross+1)*options.LineWidth*20,ModifycolorL(n,:),"filled",'DisplayName',['String_',num2str(n)]);
                end
            end
            ylabel(ax,'');
            yticks(ax,[]);
            xlabel(ax,'');
            xticks(ax,[]);
            axis(ax,'off');
        end
    end
    %% Construct Tools
    methods(Static)
        function BraidNum = BraidWord2Num(BraidWord)
            CharBraidWord = char(BraidWord);
            CharBraidWord_col = (CharBraidWord.');
            BraidNum = double(CharBraidWord_col);
            for i = 1:numel(BraidNum)
                if BraidNum(i) < 97 
                    BraidNum(i) =  BraidNum(i) - 64;
                else
                    BraidNum(i) = -BraidNum(i) + 96;
                end
                
            end
            
        end
        function Nstrings = NumBraidWord2Nstrings(DoubleBraidWord_col)
            Nstrings = max(abs(DoubleBraidWord_col)) + 1;
        end
        function [Generator,GeneratorStr] = Nstrings2Generator(Nstrings)
            TmpList =  string((1:Nstrings-1).');
            GeneratorStr = "sigma_" + TmpList;
            Generator = str2sym(GeneratorStr);
        end
        function Permutation = PermutationBraidNum(DoubleBraidWord_col)
            PermutationList = abs(DoubleBraidWord_col);
            Nstrings = max(abs(DoubleBraidWord_col)) + 1;
            Permutation(1,:) = 1:Nstrings;
            for iString = Permutation(1,:)
                Permutation(2,iString) = Braid.FinalPosition(iString,PermutationList);
            end
        end
        function oString = FinalPosition(iString,PermutationList)
            tmpString =  iString;
            for i = 1:numel(PermutationList)
                if tmpString == PermutationList(i)
                    tmpString = tmpString + 1;
                elseif tmpString == PermutationList(i) +1
                    tmpString = tmpString - 1;
                else

                end
            end
            oString = tmpString;
        end
        function [CycleDecomposition,CycleDecompositionStr,CycleLift] = Cauchy2Cycle(Permutation)
            %CyclePresentation{1} = 1;
            PermutationStore = Permutation(1,:);
            PermutationFunction = Permutation(2,:);
            Ndivide = 0;
            while ~isempty(PermutationStore)
                TriceElement = PermutationStore(1);
                Converge = false;
                RefTriceElement = TriceElement;
                Ndivide = Ndivide + 1;
                count = 0;
                while ~Converge
                    count = count + 1;
                    CycleDecomposition{Ndivide}(count) = TriceElement;
                    PermutationStore(PermutationStore == TriceElement) = [];
                    TriceElement = PermutationFunction(TriceElement);
                    Converge = TriceElement == RefTriceElement;
                end            
            end
            CycleDecompositionStr = ''; 
            for i = 1:numel(CycleDecomposition)
                CycleDecompositionStr = [CycleDecompositionStr,'('];
                CycleDecompositionStr = [CycleDecompositionStr,num2str(CycleDecomposition{i})];
                CycleDecompositionStr = [CycleDecompositionStr,')'];
                CycleLift{i} = CycleDecomposition{i};
            end
            
        end
        function Flambdak = f_lambdak(FnCell,GnCell,knCell,Ncycle,ln)
            syms lambda;
            syms k real;
            Flambdak = sym(1);
            for n = 1:Ncycle
                for j = 0:ln-1
                    Flambdak = Flambdak*(lambda-FnCell{n}(knCell{j+1})-1i*GnCell{n}(knCell{j+1}));
                end
            end

        end
    end
    methods(Static) % script
        function Permutation = PermutationBraidWord(BraidWord)
            DoubleBraidWord_col = Braid.BraidWord2Nstrings(BraidWord);
            Permutation = Braid.PermutationBraidNum(DoubleBraidWord_col);
        end
        function Nstrings = BraidWord2Nstrings(BraidWord)
            DoubleBraidWord_col = Braid.BraidWord2Nstrings(BraidWord);
            Nstrings = Braid.NumBraidWord2Nstrings(DoubleBraidWord_col);
        end
    end
end