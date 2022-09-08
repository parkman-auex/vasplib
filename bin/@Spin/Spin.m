classdef Spin < HollowKnight
    properties(Access = private,Dependent)
        s  
        ms 
    end
    properties(Hidden)
        J 
        Jz
    end
    properties
        parity = 1;
        % coe = 1;
    end
    properties
        orientation = [0 0 1];
    end
    properties(Dependent,Hidden)
        hollow;
        s2_bar;
        %parity;
        EigenVectors ;
        EigenVector;%
        %sz_bar_positive;
        %sz_bar_positive;
    end
    methods
        function SpinObj = Spin(S,Sz,coe,options)
            arguments
                S {Spin.mustBeHalfInteger(S)} = 0;
                Sz = [];
                coe = 1;
                options.orientation = [0 0 1];
                options.parity = 1;
                %options.coe = 1;
            end
            optionsCell = namedargs2cell(options);
            % call father : HollowKnight
            SpinObj = SpinObj@HollowKnight();
            %
            if isempty(S)
                SpinObj = Spin.empty([0 0]);
                return;
            end
            if isempty(Sz)
                %SpinObj.J = S;
                %SpinObj.orientation = options.orientation;
                Sz = (S:-1:-S).';
            end
            % use S = -1;Sz = nan;coe =nan; any of them repersent a empty
            % function
            
            % single
            SpinObj(1,1).J = S(1,1);
            SpinObj(1,1).orientation = options.orientation;
            if isinteger(SpinObj(1,1).J) 
                SpinObj(1,1).parity = (-1)^SpinObj(1,1).J;
            else
                SpinObj(1,1).parity = options.parity;
            end
            SpinObj(1,1).Jz = Sz(1,1);
            SpinObj(1,1).coe = coe(1,1);
            
            % multi
            [nS_row,nS_col] = size(S);
            [nSz_row,nSz_col] = size(Sz);
            [ncoe_row,ncoe_col] = size(coe);
            nSpin = max([nS_row ,nSz_row, ncoe_row] );
            nSpin_col = max([nS_col ,nSz_col, ncoe_col] );
            if nSpin == 1 && nSpin_col == 1
                return;
            else
                SpinObj = repmat(SpinObj(1,1),[nSpin,nSpin_col]);
                [S,Sz,coe] = Spin.StanderdInput(S,Sz,coe,'nrow',nSpin,'ncol',nSpin_col);
                for i = 1:numel(SpinObj)
                    SpinObj(i) = Spin(S(i),Sz(i),coe(i),optionsCell{:});
                end
            end
        end
    end
    %% overload
    methods
        function A = uminus(A)
            for i = 1:numel(A)
                A(i).coe = -A(i).coe;
            end
        end
        function C = minus(A,B,options)
           arguments
                A
                B
                options.autocontract = true;
           end
           if isa(A,'Spin') && isa(B,'Spin')
               C = A;
               C = [C,-B];
               if options.autocontract
                   C = contract(C);
               else
               end
           end
        end
        function C = plus(A,B,options)
            arguments
                A
                B
                options.autocontract = true;
            end
            if isa(A,'Spin') && isa(B,'Spin')
                C = A;
                C = [C,B];
                if options.autocontract
                    C = contract(C);
                else
                end
            end
        end
        function C = innertimes(A,B)
            [Ergebnis,S,Sz] = Spin.CGM(A.J,B.J,A.Jz,B.Jz);
            C = Spin(S,Sz,Ergebnis);
        end
        function C = eq(A,B,options)
            arguments
                A Spin
                B Spin
                options.strict = false
            end
            C = true;
            if A.s ~= B.s
                C =false;
            end
            if A.ms ~= B.ms
                C =false;
            end
            % ori?
            if ~isequal(A.orientation , B.orientation)
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
                A Spin;
                options.forgetcoe = false;
            end
            %
            B = A;%keep origin
            %cleanzero first
            A = cleanrow(A);
            % delete zero coe!!
            jmL = ([A(:).J;A(:).Jz;]).';
            coeL = [A(:).coe];
            [jmL_unique,~,~] = unique(jmL,'rows');
            if isempty(jmL_unique)
                %B = A;
                return;
            end
            if options.forgetcoe
                A = repmat(A(1),[1 length(jmL_unique)]);
                for i =1:length(jmL_unique)
                    A(i).J = jmL_unique(i,1);
                    A(i).Jz = jmL_unique(i,2);
                    A(i).coe = 1;
                end
            else
                B = repmat(A(1),[1 size(jmL_unique,1)]);
                for i = 1:size(jmL_unique,1)
                    jml_tmp = jmL_unique(i,:);
                    [~,index_ic] = ismember(jmL,jml_tmp,'rows');
                    B(1,i).J = jmL_unique(i,1);
                    B(1,i).Jz = jmL_unique(i,2);
                    B(1,i).coe = sum(coeL(logical(index_ic)));
                end
            end
            B = cleanrow(B);
            % A = Spin(nlmL_unique(nonzeroindex,2),nlmL_unique(nonzeroindex,3),coeL_unique(nonzeroindex),nlmL_unique(nonzeroindex,1));
        end
    end
    %% set 
    methods
        function SpinObj = setparity(SpinObj,paritymat)
            if isequal(size(SpinObj),size(paritymat))
                for i = 1:numel(SpinObj)
                    SpinObj(i).parity = paritymat(i);
                end
            elseif length(paritymat) == 1
                for i = 1:numel(SpinObj)
                    SpinObj(i).parity = paritymat;
                end
            end
        end
    end
    %% get
    methods
        function hollow = get.hollow(SpinObj)
            if isnan(SpinObj.coe)
                hollow = true;
            else
                hollow = false;
            end
        end
        function s = get.s(SpinObj)
            s = SpinObj.J;
        end
        function ms = get.ms(SpinObj)
            ms = SpinObj.Jz;
        end
        function s2_bar = get.s2_bar(SpinObj)
            s2_bar = SpinObj.J*(SpinObj.J+1);
        end
%         function parity = get.parity(SpinObj)
%             parity = (-1)^(SpinObj.J);
%         end
        function EigenVectors = get.EigenVectors(SpinObj)
            %for i = 1:length(SpinObj)
            i=1;
            SzM = Sz(SpinObj(i),'full',true);
            SxM = Sx(SpinObj(i),'full',true);
            SyM = Sy(SpinObj(i),'full',true);
            [W,U] = eig(SxM);
            V= W(:,(diag(U) == SpinObj(i).Jz));
            EigenVectors{i,1} = V/norm(V);
            [W,U] = eig(SyM);
            V= W(:,(diag(U) == SpinObj(i).Jz));
            EigenVectors{i,2} = V/norm(V);
            [W,U] = eig(SzM);
            V= W(:,(diag(U) == SpinObj(i).Jz));
            EigenVectors{i,3} = V/norm(V);
            %end
        end
    end
    %% rotation
    methods
        % 
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
            A_L = Spin(A.J);
            Ak = A.HollowMe;
            for i =  1:length(A_L)
                Ai = A_L(i);
                if A.J == 0
                    Ai.coe = 1;
                else
                    m1 = A.Jz;
                    m2 = Ai.Jz;
                    WignerD_single_element = (Y_l__m.d(A.J,m1,m2,beta));
                    Ai.coe = Ai.coe*exp(1i*RightorLeft*m1*alpha)*WignerD_single_element*exp(1i*RightorLeft*m2*gamma);
                    % for Y_l__m
                    % Ai.coe = conj(Ai.coe);
                end
                if immproper
                    if A.parity == 1
                        Ai.coe = Ai.coe;
                    else
                        Ai.coe = -1*Ai.coe;
                    end
                end
                if ~conjugate
                    Ai.coe = Ai.coe * A.coe;
                else
                    % test
                    Ai = TimeRerversal(Ai);
                    Ai.coe = conj(Ai.coe * A.coe);
                end
                Ak = [Ak,Ai];
            end
        end
        function U = rotation(SpinObj,rotm,rightorleft,options)
            arguments
                SpinObj Spin;
                rotm {Spin.mustBeSize(rotm,[3 3;1 5;5 1;1 4;4 1])}= diag([1 1 1]);% [alpha beta gamma]
                rightorleft = 'right';
                options.full = false;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if strcmp(rightorleft,'right')
                Coe = 1;
            else
                Coe = -1;
            end
            if isequal(size(rotm),[3 3])
                if det(rotm) == -1
                    % inversion
                    immproper = true;
                    rotm = -rotm;
                else
                    immproper = false;
                end
                %
                abc = Oper.Rotation2eul(rotm);% alpha beta gamma in ZYZ
            elseif isequal(size(rotm),[4 1]) || isequal(size(rotm),[1 4])
                if sym(rotm(end)) == -1
                    immproper = true;
                else
                    immproper = false;
                end
                if sym(abs(rotm(end))) ~= 1 
                    abc = Oper.axang2eul(rotm(1:4));
                    warning('wrong input,eular angle ZYZ right format:[alpha beta gamma det()]');
                    immproper = true;
                else
                    abc = rotm(1:3);
                end
            elseif isequal(size(rotm),[5 1]) || isequal(size(rotm),[1 5])
                if sym(abs(rotm(end))) ~= 1
                    error('wrong input,axis angle det right format:[nx ny nz theta det()]');
                end
                if sym(rotm(end)) == -1
                    immproper = true;
                else
                    immproper = false;
                end
                abc = Oper.axang2eul(rotm(1:4));
            end
            if immproper
                Invmat = ParityMat(SpinObj);
            else
                Invmat = eye(size(SpinObj,1));
            end
            if options.sym 
                abc = sym(abc);
            else
                %abc
            end
            %
            U = WignerD(SpinObj,abc,rightorleft);
            %d(j,m1,m2,seed)
            U = U *Invmat;
            if isa(U,'sym')
                U = simplify(U);
            end
        end
        function Trmat = Tr(SpinObj)
            Trmat = InnerProduct(SpinObj,TimeRerversal(SpinObj));
        end
        function Invmat = ParityMat(SpinObj)
            Invmat = InnerProduct(SpinObj,Inversion(SpinObj));
        end
        function WigerDmat = WignerD(SpinObj,abc,rightorleft)
            arguments
                SpinObj Spin;
                abc;
                rightorleft = 'right';
            end
            nSpin = size(SpinObj,1)   ;
            WigerDmat =sym( zeros(nSpin));
            for i = 1:nSpin
                for j = 1:nSpin
                    Matelement = WignerD_single(SpinObj(i,:),SpinObj(j,:),abc,rightorleft);
                    if Matelement~=sym(0)
                        WigerDmat(i,j) = Matelement;
                    end
                end
            end
        end
        function Matelement = WignerD_single(SpinObj1,SpinObj2,abc,rightorleft)
            arguments
                SpinObj1 Spin;
                SpinObj2 Spin;
                abc;
                rightorleft = 'right';
            end
            assert(isrow(SpinObj1));
            assert(isrow(SpinObj2));
            alpha = sym(abc(1));
            beta = sym(abc(2));
            gamma = sym(abc(3));
            Matelement = sym(0);
            if strcmp(rightorleft,'right')
                for i = 1:length(SpinObj1)
                    for j = 1:length(SpinObj2)
                        if SpinObj1.J ~= SpinObj2.J
                            continue;
                        end
                        m1 = SpinObj1.Jz;
                        m2 = SpinObj2.Jz;
                        WignerD_single_element = sym(Y_l__m.d(SpinObj1.J,m1,m2,beta));
                        Matelement = Matelement + exp(1i*m1*alpha)*WignerD_single_element*exp(1i*m2*gamma);
                    end
                end
            end
        end
        function U = rotation2(SpinObj,rotm,rightorleft,options)
            % not recommand delete later!!
            arguments
                SpinObj Spin;
                rotm {Spin.mustBeSize(rotm,[3 3;1 5;5 1;1 4;4 1])} = diag([1 1 1]);% [alpha beta gamma]
                rightorleft = 'right';
                options.full = false;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if strcmp(rightorleft,'right')
                Coe = 1;
            else
                Coe = -1;
            end
            if isequal(size(rotm),[3 3])
                if det(rotm) == -1
                    % inversion
                    immproper = true;
                    rotm = -rotm;
                else
                    immproper = false;
                end
                %
                abc = Oper.Rotation2eul(rotm);% alpha beta gamma in ZYZ
            elseif isequal(size(rotm),[4 1]) || isequal(size(rotm),[1 4])
                if sym(rotm(end)) == -1
                    immproper = true;
                else
                    immproper = false;
                end
                if sym(abs(rotm(end))) ~= 1 
                    abc = Oper.axang2eul(rotm(1:4));
                    warning('wrong input,eular angle ZYZ right format:[alpha beta gamma det()]');
                    immproper = true;
                else
                    abc = rotm(1:3);
                end
            elseif isequal(size(rotm),[5 1]) || isequal(size(rotm),[1 5])
                if sym(abs(rotm(end))) ~= 1
                    error('wrong input,axis angle det right format:[nx ny nz theta det()]');
                end
                if sym(rotm(end)) == -1
                    immproper = true;
                else
                    immproper = false;
                end
                abc = Oper.axang2eul(rotm(1:4));
            end
            if immproper
                Invmat = ParityMat(SpinObj);
            else
                Invmat = eye(size(SpinObj,1));
            end
            if options.sym 
                abc = sym(abc);
            else
                %abc
            end
            %
            alpha = abc(1);
            beta = abc(2);
            gamma = abc(3);
            U = expm(Coe*1i*(...
                alpha*Sz(SpinObj,optionsCell{:})+...
                beta*Sy(SpinObj,optionsCell{:})+...
                gamma*Sz(SpinObj,optionsCell{:}) ...
                ));
            U = U *Invmat;
            if isa(U,'sym')
                U = simplify(U);
            end
        end
    end
    %% operator
    methods
        function SpinObj = CG(Spinobj1,Spinobj2,options)
            
        end
        function SpinObj = TimeRerversal(SpinObj)
            % \Theta|j m\rangle=(-1)^{j-m}|j,-m\rangle
            for i = 1:numel(SpinObj)
                SpinObj(i).coe = (-1)^(SpinObj(i).J-SpinObj(i).Jz)*SpinObj(i).coe;
                SpinObj(i).Jz = -SpinObj(i).Jz;
            end
        end
        function SpinObj = Inversion(SpinObj)
            for i =1:numel(SpinObj)
                SpinObj(i).coe = SpinObj(i).coe * SpinObj(i).parity;
            end
        end
        function SpinObj = SzOper(Spinobj)
            SpinObj = Spinobj;
            for i = 1:numel(Spinobj)
                SpinObj(i).coe = SpinObj(i).coe * SpinObj(i).coe*Spinobj(i).Jz;
            end
            SpinObj = SpinObj.contract();
        end
        function SpinObj = SxOper(Spinobj)
            SpinObj = 1/2.*(SplusOper(Spinobj)+SminusOper(Spinobj));
            SpinObj = SpinObj.contract();
        end
        function SpinObj = SyOper(Spinobj)
            SpinObj = 1/2i.*(SplusOper(Spinobj) - SminusOper(Spinobj));
            SpinObj = SpinObj.contract();
        end
        function SpinObj = SplusOper(Spinobj)
            SpinObj = Spinobj;
            for i = 1:numel(Spinobj)
                SpinObj(i).coe = SpinObj(i).coe * sqrt(Spinobj(i).s2_bar-Spinobj(i).Jz*(Spinobj(i).Jz+1));
                SpinObj(i).Jz = Spinobj(i).Jz+1;
            end
            SpinObj = SpinObj.contract();
        end
        function SpinObj = SminusOper(Spinobj)
            SpinObj = Spinobj;
            for i = 1:numel(Spinobj)
                SpinObj(i).coe = SpinObj(i).coe * sqrt(Spinobj(i).s2_bar-Spinobj(i).Jz*(Spinobj(i).Jz-1));
                SpinObj(i).Jz = Spinobj(i).Jz-1;
            end
            SpinObj = SpinObj.contract();
        end
    end
    %% operator mat
    methods
        %name as S
        function SzM = Sz(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = false;
                options.sym = true;
            end
            if options.full
                SpinObj = SpinObj.BasisGen();
            end
            optionsCell = namedargs2cell(options);
            nSpin = length(SpinObj);
            if options.full
                if nSpin ~= SpinObj(1).J*2+1
                    error('!');
                end
            end
            if options.sym
                S = sym(SpinObj(1).J);
                SzM_L = repmat(S,[nSpin-1 1]);
            else
                S = SpinObj.J;
                SzM_L = repmat(S,[nSpin-1 1]);
            end
            if options.full
                for a =1:nSpin
                    %b = a+1;
                    SzM_L(a) =2*(S+1-a)/2;
                end
                SzM = diag(SzM_L);
            else %depend on basis
                SzM = InnerProduct(SpinObj,SzOper(SpinObj),'sym',options.sym );
            end
        end
        function SyM = Sy(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = false;
                options.sym = true;
            end
            if options.full
                SpinObj = SpinObj.BasisGen();
            end
            optionsCell = namedargs2cell(options);
            nSpin = length(SpinObj);
            if options.full
                if nSpin ~= SpinObj(1).J*2+1
                    error('!');
                end
            end
            if options.sym
                S = sym(SpinObj(1).J);
                SyM_L = repmat(S,[nSpin-1 1]);
            else
                S = SpinObj.J;
                SyM_L = repmat(S,[nSpin-1 1]);
            end
            if options.full
                for a =1:nSpin-1
                    b = a+1;
                    SyM_L(a) =-1i*sqrt((S+1)*(a+b-1)-a*b)/2;
                end
                SyM = diag(SyM_L,1)+diag(conj(SyM_L),-1);
            else %depend on basis
                SyM = 1/2i*(Splus(SpinObj,optionsCell{:})-Sminus(SpinObj,optionsCell{:}));
            end
        end
        function SxM = Sx(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = false;
                options.sym = true;
            end
            if options.full
                SpinObj = SpinObj.BasisGen();
            end
            optionsCell = namedargs2cell(options);
            nSpin = length(SpinObj);
            if options.full
                if nSpin ~= SpinObj(1).J*2+1
                    error('!');
                end
            end
            if options.sym
                S = sym(SpinObj(1).J);
                SxM_L = repmat(S,[nSpin-1 1]);
            else
                S = SpinObj.J;
                SxM_L = repmat(S,[nSpin-1 1]);
            end
            if options.full
                for a =1:nSpin-1
                    b = a+1;
                    SxM_L(a) =sqrt((S+1)*(a+b-1)-a*b)/2;
                end
                SxM = diag(SxM_L,1)+diag(conj(SxM_L),-1);
            else %depend on basis
                SxM = 1/2*(Splus(SpinObj,optionsCell{:})+Sminus(SpinObj,optionsCell{:}));
            end
        end
        %name as L
        function LzM = Lz(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            LzM = Sz(SpinObj,optionsCell{:});
        end
        function LxM = Lx(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            LxM = Sx(SpinObj,optionsCell{:});          
        end
        function LyM = Ly(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            LyM = Sy(SpinObj,optionsCell{:});
        end
        %name as J
        function JzM = JZ(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            JzM = Sz(SpinObj,optionsCell{:});
        end
        function JxM = JX(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            JxM = Sx(SpinObj,optionsCell{:});
        end
        function JyM = JY(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            JyM = Sy(SpinObj,optionsCell{:});
        end
        % Ladder function
        function Splus = Splus(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if options.full
                Splus = Sx(SpinObj,optionsCell{:})+1i*Sy(SpinObj,optionsCell{:});
            else
                Splus = InnerProduct(SpinObj,SplusOper(SpinObj),'sym',options.sym );
            end
        end
        function Sminus = Sminus(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if options.full
                Sminus = Sx(SpinObj,optionsCell{:})-1i*Sy(SpinObj,optionsCell{:});
            else
                Sminus = InnerProduct(SpinObj,SminusOper(SpinObj),'sym',options.sym );
            end
        end
        %name as L
        function Lplus = Lplus(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            Lplus = Lx(SpinObj,optionsCell{:})+1i*Ly(SpinObj,optionsCell{:});
        end
        function Lminus = Lminus(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            Lminus = Lx(SpinObj,optionsCell{:})-1i*Ly(SpinObj,optionsCell{:});
        end
        %name as J
        function Jplus = Jplus(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            Jplus = JX(SpinObj,optionsCell{:})+1i*JY(SpinObj,optionsCell{:});
        end
        function Jminus = Jminus(SpinObj,options)
            arguments
                SpinObj Spin;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            Jminus = JX(SpinObj,optionsCell{:})-1i*JY(SpinObj,optionsCell{:});
        end
    end

    %% Modify
    methods
        function SpinObj = BasisGen(Spinobj,force)
            if nargin <2
                force = true;
            end
            if force
                count = 0;
                for Sz = Spinobj(1).J:-1:-Spinobj(1).J
                    count = count + 1;
                    SpinObj(count) = Spin(Spinobj(1).J,Sz,'orientation',Spinobj(1).orientation);
                end
            end
        end
        function SpinObj = PerM(Spinobj,permutation)

        end
    end
    %% pretty
    methods
        function disp(SpinObj)
            strlist = string(SpinObj);
            for i = 1:size(strlist,1)
                fprintf(strlist(i,1));
            end
        end
        function Output = string(SpinObj)
            Output = pretty(SpinObj);
        end
        function Output = pretty(SpinObj,options)
            arguments
                SpinObj Spin
                options.output char{mustBeMember(options.output,{'str','sym','latex'})} = 'str';
                options.formatspec = '%5.3f';
            end
            if isempty(SpinObj)
                Output = "";
                return;
            end
            %
            optionsCell = namedargs2cell(options);
            rightbraket = '\x27E9';
            %
            if length(SpinObj) == 1
                if SpinObj.hollow
                    Output = "";
                    return;
                end
                switch options.output
                    case 'str'
                        if SpinObj.coe == 1 || SpinObj.coe == sym(1)
                            strcoe = "";
                        elseif SpinObj.coe == -1 || SpinObj.coe == sym(-1)
                            strcoe = "-";
                        elseif ~isreal(SpinObj.coe) 
                            strcoe = string(sym(SpinObj.coe));
                        elseif isa(SpinObj.coe,'sym')
                            strcoe = string(sym(SpinObj.coe));
                        else
                            strcoe = num2str((SpinObj.coe),options.formatspec);
                        end
                        Output = strcoe +"|"+string(sym(SpinObj.J))+","+string(sym(SpinObj.Jz))+rightbraket+" ";
                    case 'latex'
                    case 'sym'
                end
            elseif length(SpinObj) > 1
                switch options.output
                    case 'str'
                        for i = 1:size(SpinObj,1)
                            Output(i,1) = pretty(SpinObj(i,1),optionsCell{:});
                            for j =2:size(SpinObj,2)
                                Output(i,1) = Output(i,1) +"+ " +  pretty(SpinObj(i,j),optionsCell{:});
                            end
                            Output(i,1) = Output(i,1)+"\n";
                        end
                    case 'latex'
                    case 'sym'
                end


            else
            end
        end
    end

    %% math
    methods(Static)
        function [Ergebnis,j3,m3] = CGM(j1,j2,m1,m2)
            % Evaluation of all Clebsch-Gordan coefficients for fixed M and
            % m1 m2 j1 j2
            % call: CB =ClebschGordan(j1,j2,m3)
            % Example: j1 = 1/2 j2=1 m3 = 1/2 m1 =1/2 m2 =1
            % based on Vers. 2.2 18.08.2016 ClebschGordan.m
            % Copyright Wolfgang Schweizer
            % Simulation physikalischer Systeme - Computational Physics mit MATLAB
            % 2016
            % Vers. 1.1
            % 18.10.2019
            if j1 ==-1 ||j2 == -1|| isnan(m1) || isnan(m2)
                Ergebnis = nan;j3 = -1;m3 = nan;
                return;
            end
            m3 = m1+m2;
            j3 = abs(j1-j2):j1+j2;
            j3 = j3(j3>=abs(m3));   % allowed j3 values
            for n=1:length(j3)
                Ergebnis(1,n) = Spin.CGsingle(j1,j2,j3(n),m1,m2,m3);
            end
        end
        function Ergebnis = CGsingle(j1,j2,j3,m1,m2,m3)
            % Evaluation of Clebsch-Gordan coefficient based on
            %               all angular values
            % CG = CGsingle(j1,j2,j3,m1,m2,m3)

            % based on Vers. 2.2 18.08.2016
            % Wolfgang Schweizer
            % Simulation physikalischer Systeme - Computational Physics mit MATLAB
            % Vers. 1.1 18.10.21019

            Ergebnis = [];
            %
            vorfakj3 = [j1+j2-j3,j3+j1-j2,j3+j2-j1];
            vorfakj3 = factorial(vorfakj3);
            vorfakj3 = sqrt((2*j3+1)*prod(vorfakj3)/factorial(j3+j1+j2+1));
            sums = [];
            % test via Regge symbol
            % Racah 用代数方法得出了克莱布希－高登系数的有限级数表达式[3]。
            % Giulio Racah. Theory of Complex Spectra. II. Phys. Rev.: 438. doi:10.1103/PhysRev.62.438.
            regge = [-j1+j2+j3   , j1-j2+j3   , j1+j2-j3    ; ...
                j1 - m1    ,   j2-m2    ,   j3+m3  ; ...
                j1 + m1    ,   j2+m2    ,   j3-m3 ];
            if sum(regge(:) < 0) || sum(sum(regge)~=j1+j2+j3) || sum(sum(regge,2)~=j1+j2+j3)
                disp('Null!! Incorrect values')
                Ergebnis = [Ergebnis;0];
                %return
            else
                %[j1,m1;j2,m2;j3,m3];
                lauf1 = [j1+j2-j3,j1-m1,j2+m2,j3-j2+m1,j3-j1-m2];
                s=abs(min(lauf1));
                laufwhile = logical(1);
                notstop=1;
                vorfak = [j1+m1,j1-m1,j2+m2,j2-m2,j3+m3,j3-m3];
                %
                if sum(vorfak<0)
                    disp('Erroneous Prefactor')
                else
                    vorfak = factorial(vorfak);
                    %vorfak = sqrt(prod(vorfak));
                    vorfak = exp(sum(log(vorfak)*0.5));
                end
                %
                while laufwhile
                    test = [j1+j2-j3-s,j1-m1-s,j2+m2-s,j3-j2+m1+s,j3-j1-m2+s];
                    testok=min(test); %#ok<*NASGU>
                    zwi1 = [s,test];
                    zwi1 = factorial(zwi1);
                    %zwi1 = prod(zwi1);
                    zwi1 = exp(sum(log(zwi1)));
                    sums(notstop) = (-1)^s * vorfak / zwi1;  %#ok<*AGROW>
                    s=s+1;
                    test = [j1+j2-j3-s,j1-m1-s,j2+m2-s,j3-j2+m1+s,j3-j1-m2+s];
                    testok=min(test);
                    notstop = notstop+1;
                    if notstop > 100
                        disp('while-loop run 100 times - error?')
                        return
                    end
                    laufwhile = testok>=0;
                end
                Summe = sum(sums);
                Ergebnis = [Ergebnis;vorfakj3*sym(Summe)];
                %ClebschGordanKoeff = array2table(Ergebnis);
                %ClebschGordanKoeff.Properties.VariableNames ={'j1' 'j2' 'J' 'm1' 'm2' 'M' 'CG'};
            end

        end
    end
    %% tool
    methods(Static)
        function mustBeHalfInteger(S)
            mustBeInteger(S*2);
        end
        function mustBeSize(a,b)
            % Test for size
            if size(b,1) ==1
                if sum(~(size(a)==b))
                    eid = 'Size:notRequired';
                    msg = 'Inputs must have size of '+mat2str(b);
                    throwAsCaller(MException(eid,msg))
                end
            elseif   size(b,1) ==2
                if sum(~(size(a)==b(1,:))) && sum(~(size(a)==b(2,:)))
                    eid = 'Size:notRequired';
                    msg = 'Inputs must have size of '+mat2str(b);
                    throwAsCaller(MException(eid,msg))
                end
            else
            end
        end
        %         function T = transferoptions(options)
        %             T_tmp = namedargs2cell(options);
        %             T = T_tmp{:};
        %         end
    end
end