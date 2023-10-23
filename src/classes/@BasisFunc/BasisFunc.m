classdef BasisFunc < HollowKnight
    properties
        BForb = [];
        BFnum =[];
        BFuncL = [];
        spin = [];
    end
    properties(Hidden,Dependent)
        hollow;
    end
    properties(Access = private,Dependent)
        FuncNum ;
        BFclassL ;
        % Expression;
    end
     
    properties(Hidden)
        parity;
        Rm = [1 0 0;0 1 0;0 0 1];
    end
    %% construction inner
    methods
        function BasisFunction = BasisFunc(BFuncLOrigin,spin,BFnum,coe,orb,options)
            % BasisFunction = BasisFunc({Spin(1/2)})
            % BasisFunction = BasisFunc(Spin(1/2))
            arguments
                BFuncLOrigin = Qnum().HollowMe;
                spin  Spin = Spin(0,0) ;
                BFnum = [];
                coe = sym([]);
                orb = [];
                options.raw = true
                options.fast = true;
            end
            optionsCell = namedargs2cell(options);
            if isempty(BFuncLOrigin)
                BasisFunction = BasisFunc.empty([1 0]);
                return;
            end
            if isa(BFuncLOrigin,'double') 
                SizeBasisFunc = BFuncLOrigin;
                if isempty(BFuncLOrigin)
                    BasisFunction = BasisFunc.empty([1 0]);
                elseif length(SizeBasisFunc) == 1
                    BasisFunction = BasisFunc([SizeBasisFunc 1]);
                    return;
                else
                    BasisFunction1 = BasisFunc();
                    BasisFunction = repmat(BasisFunction1,SizeBasisFunc);
                    return;
                end
            elseif isa(BFuncLOrigin,'vasplib')
                BasisFunction = BasisFunc.BasisFunction(BFuncLOrigin);
                return;
            end
            if isa(BFuncLOrigin,'Spin') ||isa(BFuncLOrigin,'jm')
                spin = BFuncLOrigin;
                BFuncLOrigin = sym(1);
            end
            if isempty(BFnum)                % coe
                BFnum = 1;
            end
            if isempty(coe)
                coe = sym(1);
            end
            BasisFunction(1,1).BFuncL = BFuncLOrigin(1,1);
            BasisFunction(1,1).spin = spin(1,1);
            BasisFunction(1,1).BFnum = BFnum(1,1);
            BasisFunction(1,1).coe = coe(1,1);
            if ~isempty(orb)
                BasisFunction(1,1).BForb = orb(1,:);
            end
            if isscalar(BFuncLOrigin) && isscalar(spin)
                return;
            end
            sizeBFuncL = size(BFuncLOrigin);
            sizespin = size(spin);
            if isequal(sizeBFuncL,sizespin)
            elseif isvector(spin) && isvector(BFuncLOrigin)
                spin = repmat(spin,sizeBFuncL);
                BFuncLOrigin = repmat(BFuncLOrigin,sizespin);
            else
                try
                    if sizeBFuncL(2) > sizespin(1)
                        spin = repmat(spin,sizeBFuncL./sizespin);
                    else
                        BFuncLOrigin = repmat(BFuncLOrigin,sizespin./sizeBFuncL);
                    end
                catch
                    error('wrong input!')
                end
            end
            sizeBF = size(BFuncLOrigin);
            if ~isequal(sizeBF,size(BFnum))
                BFnum = repmat(BFnum,sizeBF./size(BFnum));
            end
            if ~isequal(sizeBF,size(coe))
                coe = repmat(coe,sizeBF./size(coe));
            end
            BasisFunction = repmat(BasisFunction(1,1),sizeBF);
            for i = 1:numel(BasisFunction)
                BasisFunction(i) = BasisFunc(BFuncLOrigin(i),spin(i),BFnum(i),coe(i),optionsCell{:});
            end
            if ~isempty(orb)
                norbL = size(orb,1);
                if norbL == 1
                    for i = 1:numel(BasisFunction)
                        BasisFunction(i).BForb = orb(1,:);
                    end
                elseif norbL == numel(BasisFunction)
                    for i = 1:numel(BasisFunction)
                        BasisFunction(i).BForb = orb(i,:);
                    end
                elseif norbL == size(BasisFunction,1)
                    for i = 1:size(BasisFunction,1)
                        for j = 1:size(BasisFunction,1)
                            BasisFunction(i,j).BForb = orb(i,:);
                        end
                    end

                else

                end
            end
        end
    end
    %% construction outer
    methods(Static)
        function BasisFunction = BasisFunction(vasplibobj)
            switch class(vasplibobj)
                case {'vasplib','HR','Htrig','HK'}
                    % Qnum
                    [BFuncLOrigin,S,SzL] = Qnum.QnumL(vasplibobj);
                    spinL = Spin(S,SzL);
                    BasisFunction = BasisFunc(BFuncLOrigin,spinL,1,1,vasplibobj.orbL);
                    for i = 1:numel(BasisFunction)
                        BasisFunction(i).Rm = vasplibobj.Rm;
                    end
            end

        end
    end
    %% get
    methods
        function hollow = get.hollow(BasisFunction)
            try
                hollowBFuncL = BasisFunction.BFuncL.hollow;
            catch
                hollowBFuncL = false;
            end
            if hollowBFuncL || BasisFunction.spin.hollow || isnan(BasisFunction.coe)
                hollow = true;
            else
                hollow = false;
            end
        end
        function FuncNum = get.FuncNum(BasisFunction)
            FuncNum = length(BasisFunction.BFuncL);
        end
        function BFclassL = get.BFclassL(BasisFunction)
            NumFunc = BasisFunction.FuncNum;
            BFclassL(NumFunc) = string(class(BasisFunction.BFuncL{NumFunc}));
            for i = 1:BasisFunction.FuncNum-1
                BFclassL(i) = string(class(BasisFunction.BFuncL{i}));
            end
        end
%         function Expression = get.Expression(BasisFunction)
%             Expression = string(BasisFunction);
%         end
    end
    %% disp
    methods
        function dispAll(BasisFunction)
            for i = 1:size(BasisFunction,1)
                fprintf('============================================\n')
                %fprintf('The %d / %d th BasisFunction:\n',i,numel(BasisFunction));
                fprintf(string(BasisFunction(i,1).coe));
                fprintf(string(BasisFunction(i,1).BFuncL));
                fprintf(string(BasisFunction(i,1).spin));
                for j = 2:size(BasisFunction,2)
                    fprintf(' + ');
                    fprintf(string(BasisFunction(i,j).coe));
                    fprintf(string(BasisFunction(i,j).BFuncL));
                    fprintf(string(BasisFunction(i,j).spin));
                end
                fprintf('\n')
            end
            fprintf('============================================\n')
        end
        function disporbL(BasisFunction)
            disp(reshape([BasisFunction(:,1).BForb],[],3));
        end
    end
    %% overload
    methods
        function C = eq(A,B,options)
            arguments
                A BasisFunc;
                B BasisFunc;
                options.strict = false;
                options.spin = false;
                options.orb = true;
                options.BFuncL =true;
            end
            C = true;
            if class(A.BFuncL) ~= class(B.BFuncL)
                C =false;
                return;
            end
            if options.spin 
                if ~eq(A.spin,B.spin,'strict',options.strict)
                    C =false;
                    return;
                end
            end
            if options.orb
                if isempty(A.BForb) &&  ~isempty(B.BForb) 
                    C =false;
                    return;
                end
                if ~isempty(A.BForb) &&  isempty(B.BForb)
                    C =false;
                    return;
                end
                if ~Oper.allclose(A.BForb,B.BForb)
                    C =false;
                    return;
                end
            end
            if options.BFuncL
                if ~eq(A.BFuncL,B.BFuncL,'strict',options.strict)
                    C =false;
                    return;
                end
            end
        end
        function C = innertimes(A,B,options)
            arguments
                A BasisFunc;
                B BasisFunc;
                options.forgetcoe = false;
                options.fast = true;
                options.hybird = false;
                options.spincoupled = false;
                options.orbcoupled = false;
                options.raw = true
            end
            optionsCell = namedargs2cell(options);
            if ~options.hybird && ~options.spincoupled && ~options.orbcoupled
                coeLtmp0 = A.coe*B.coe;
                BFuncLtmp = (A.BFuncL .* B.BFuncL);
                [coeLtmp,BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp);
                spinL = A(1).spin;
                orbL = A(1).BForb;
                C = BasisFunc(BFuncLtmp,spinL,1,coeLtmp0*coeLtmp,'orbL',orbL);
            end
        end
    end
    %% contract
    methods
        % we can overload contract too
        function B = contractrow(A,options)
            arguments
                A BasisFunc;
                options.forgetcoe = false;
                options.fast = true;
                options.hybird = false;
                options.spincoupled = false;
                options.orbcoupled = false;
                options.raw = true;
                options.sym = false;
                options.conjugate = false;
                options.antisymmetry = false;
            end
            % 
            B = A;%keep origin
            % for hybird 
            if options.hybird
                error('not support yet');
                return
            end
            % for cell  
            if iscell(A(1).BFuncL)
                error('not support yet');
                return;
            end
            A = cleanrow(A);
            if  options.spincoupled
            end
            if options.orbcoupled 
            end
            if ~options.spincoupled && ~options.orbcoupled
                BFuncLtmp = ([A.BFuncL]);
                coeLtmp = ([A.coe]); % num?
                % giveback coe
                BFuncLtmp = BasisFunc.introduce_coe(BFuncLtmp,coeLtmp);
                % contract row
                BFuncLtmp = contractrow(BFuncLtmp);
                % extrack coe
                [coeLtmp,BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp);
                switch class(BFuncLtmp)
                    case 'Qnum'
                        spinL = A(1).spin;
                    otherwise
                        spinL = A(1).spin;
                end
                % if spin coupled?
                orbL = A(1).BForb;
                B = BasisFunc(BFuncLtmp,spinL,1,coeLtmp,orbL);
            end
            B = cleanrow(B);
        end
    end
    methods(Static)
        function [BFuncL] = introduce_coe(BFuncL,coeL)
            for i =1:numel(BFuncL)
                BFuncL(i).coe = BFuncL(i).coe*coeL(i);
            end
        end
        function [coeL,BFuncL] = extract_coe(BFuncL,options)
            arguments
                BFuncL
                options.sym = false;
                options.vpalevel = 6;
            end
            coeL = ones(size(BFuncL),class(BFuncL(1).coe));
            oneterm  = ones(1,1,class(BFuncL(1).coe));
            for i =1:numel(BFuncL)
                coeL(i) = BFuncL(i).coe;
                BFuncL(i).coe =  oneterm;
            end
            % Accuracy?
            if options.sym
                
            else
                coeL = round(coeL,options.vpalevel);
            end
        end
    end
    %% rotation
    methods
        function U = rotation(A,Rc,Rf,tf,optionsConvection,optionsOper,optionsRm,options)
            % Input arguments
            %   A       Basisfunction      (bandlimit+1)^2 x 1 complex array, spherical
            %                           harmonics coefficients.
            %   rotm    double      4 x 1 array or 3 x 3 array, rotation in quaternion
            %                           representation or in rotation matrix
            %                           representation, respectively.
            %
            arguments
                A BasisFunc;
                Rc {Spin.mustBeSize(Rc,[3 3])}= diag([1 1 1]);%
                Rf {Spin.mustBeSize(Rf,[3 3])}= diag([1 1 1]);%
                tf {Spin.mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);%
                optionsConvection.rightorleft {mustBeMember(optionsConvection.rightorleft,{'right','left'})}= 'right';
                optionsOper.Oper = [];
                optionsRm.Rm = POSCAR_read;
                options.sym = false;
                options.conjugate = false;
                options.antisymmetry = false;
                options.forgetcoe = false;
                options.fast = true;
                options.hybird = false;
                options.spincoupled = false;
                options.orbcoupled = false;
                options.raw = true;
                options.vpalevel = 6;
                options.center = [0,0,0];
            end
            rightorleft = optionsConvection.rightorleft ;
            optionsCell = namedargs2cell(options);
            if ~isempty(optionsOper.Oper)
                Rf = optionsOper.Oper.Rf;
                if isempty(Rf)
                    optionsOper.Oper =optionsOper.Oper.attachRm(optionsRm.Rm);  
                end
                Rf = optionsOper.Oper.Rf;
                Rc = optionsOper.Oper.R;
                tf = optionsOper.Oper.tf;
                options.conjugate = optionsOper.Oper.conjugate;
                options.antisymmetry = optionsOper.Oper.antisymmetry;
                optionsCell = namedargs2cell(options);
                Am = rotate(A,Rc,Rf,tf,rightorleft,optionsCell{:});
            else
                Am = rotate(A,Rc,Rf,tf,rightorleft,'Oper',optionsOper.Oper,optionsCell{:});
            end
            U  = InnerProduct(Am,A,'sym',options.sym);
        end
        function Am = rotate(A,Rc,Rf,tf,rightorleft,optionsOper,options)
            % Input arguments
            %   A         Basisfunction      (bandlimit+1)^2 x 1 complex array, spherical
            %                           harmonics coefficients.
            %   rotm    double      4 x 1 array or 3 x 3 array, rotation in quaternion
            %                           representation or in rotation matrix
            %                           representation, respectively.
            %
            arguments
                A BasisFunc;
                Rc {Spin.mustBeSize(Rc,[3 3])}= diag([1 1 1]);%
                Rf {Spin.mustBeSize(Rf,[3 3])}= diag([1 1 1]);%
                tf {Spin.mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);%
                rightorleft = 'right';
                optionsOper.Oper = [];
                options.sym = false;
                options.conjugate = false;
                options.antisymmetry = false;
                options.forgetcoe = false;
                options.fast = true;
                options.hybird = false;
                options.spincoupled = false;
                options.orbcoupled = false;
                options.raw = true;
                options.vpalevel = 6;
                options.center = [0,0,0];
            end
            if ~isempty(optionsOper.Oper)
                Rc = optionsOper.Oper.R;
                Rf = optionsOper.Oper.Rf;
                tf = optionsOper.Oper.tf;
                options.conjugate = optionsOper.Oper.conjugate;
                options.antisymmetry = optionsOper.Oper.antisymmetry;
                optionsCell = namedargs2cell(options);
                Am = rotate(A,Rc,Rf,tf,rightorleft,optionsCell{:});
                return
            end
            optionsCell = namedargs2cell(options);
            Am = rotaterow(A(1,:),Rc,Rf,tf,rightorleft,optionsCell{:});
            for i = 2:size(A,1)
                Am = [Am;rotaterow(A(i,:),Rc,Rf,tf,rightorleft,optionsCell{:})];
            end
        end
        function A_Lj = rotaterow(A,Rc,Rf,tf,rightorleft,options)
            arguments
                A BasisFunc;
                Rc {Spin.mustBeSize(Rc,[3 3])}= diag([1 1 1]);%
                Rf {Spin.mustBeSize(Rf,[3 3])}= diag([1 1 1]);%
                tf {Spin.mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);%
                rightorleft = 'right';
                options.sym = true;
                options.conjugate = false;
                options.antisymmetry = false;
                options.forgetcoe = false;
                options.fast = true;
                options.hybird = false;
                options.spincoupled = false;
                options.orbcoupled = false;
                options.raw = true;
                options.vpalevel = 6;
                options.center = [0,0,0];
            end
            optionsCell = namedargs2cell(options);
            % for hybird
            if options.hybird
                % the classtype of BFuncL should be same;
                error('not support yet');
                return
            else
                BFuncLtmp = ([A.BFuncL]);
                coeLtmp = ([A.coe]); % num?
                % giveback coe
                BFuncLtmp = BasisFunc.introduce_coe(BFuncLtmp,coeLtmp);
            end
            % for cell
            if iscell(A(1).BFuncL)
                error('not support yet');
                return;
            end
            if ~options.orbcoupled
                orbL = BasisFunc.rotation_orb(A(1).BForb,Rf.',tf,optionsCell{:});
                if ~options.spincoupled
                    if isa(A(1).BFuncL,'Qnum')
                        BFuncLtmp = rotaterow(BFuncLtmp,Rc,rightorleft,...
                            'sym',options.sym,'antisymmetry',options.antisymmetry,'conjugate',options.conjugate);
                        % contract row
                        BFuncLtmp = contractrow(BFuncLtmp);
                        %
                        spinL = Spin([BFuncLtmp.s],[BFuncLtmp.sz]);
                        % extract coe
                        [coeLtmp,BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp,'sym',options.sym,'vpalevel',options.vpalevel);
                        A_Lj = BasisFunc(BFuncLtmp,spinL,1,coeLtmp,orbL);
                    else
                        BFuncLtmp = rotaterow(BFuncLtmp,Rc,rightorleft,...
                            'sym',options.sym,'antisymmetry',options.antisymmetry,'conjugate',options.conjugate);
                        spinL = rotaterow([A.spin],Rc,rightorleft,...
                            'sym',options.sym,'antisymmetry',options.antisymmetry,'conjugate',options.conjugate);
                        % contract row
                        BFuncLtmp = contractrow(BFuncLtmp);
                        spinL = contractrow(spinL);
                    end
                else

                end
            else
                % not imply yet
                error('not imply yet');
            end
        end
        function BasisFunction = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
            switch class(A.BFuncL)
                case 'QnumL'

                otherwise
                    BFuncL = rotateinner(A.BFuncL,abc,RightorLeft,immproper,conjugate,antisymmetry);
                    SpinL = rotateinner(A.spin,abc,RightorLeft,immproper,conjugate,antisymmetry);
                    BForb = BasisFunc.rotation_orb(A.BForb,R,t,options);
                    BasisFunction = BasisFunc();
            end
        end
    end
    methods(Static)
        function [BFuncL,coeL] = rotation_func(BFuncL,R,t,options)
            arguments
                BFuncL          
                R
                t
                options.conjugate logical = false ;
                options.Rm  = [1 0 0;0 1 0;0 0 1] ;
            end
            optionsCell = namedargs2cell(options);
            for i =numel(BFuncL)
                if isa(BFuncL{i},'sym')
                    
                else
                    BFuncL{i} = rotate(BFuncL{i},R,t,optionsCell{:});
                end
            end
        end
        function spinL = rotation_spin(spin,R,t,options)

        end
        function BForb = rotation_orb(BForb,Rf,tf,options)
            arguments
                BForb 
                Rf {Spin.mustBeSize(Rf,[3 3])}= diag([1 1 1]);% 
                tf {Spin.mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);% 
                options.sym = true;
                options.conjugate = false;
                options.antisymmetry = false;
                options.forgetcoe = false;
                options.fast = true;
                options.hybird = false;
                options.spincoupled = false;
                options.orbcoupled = false;
                options.raw = true;
                options.vpalevel = 6;
                options.center = [0,0,0];
            end
            if isempty(BForb)
                return;
            end
            BForb = (BForb-options.center)*Rf.'+options.center + tf;
            % to homecell
            for i = 1:3
               BForb(i) =mod(BForb(i),1);
            end
        end
    end
    methods(Static)

    end
end