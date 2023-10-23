classdef HollowKnight %(Abstract) Reconstuct HollowKnight and group as Abstract class! 
    % There are a bunch of class need hollow obj 
    % take it as a father class
    % SubClass of HollowKnight has properties: coe and hollow
    % coe = nan means this obj is only a hollow obj
    properties
        coe {mustBeA(coe,["double","sym"])}= 1
    end
    properties(Abstract,Hidden)
        hollow;
    end
    methods
        function obj = HollowKnight()
            % empty 
            % for convenient construct subclass
        end
    end
    %% get
    methods
%         function hollow = get.hollow(HollowKnightObj)
%             if isnan(HollowKnightObj.coe)
%                 hollow = true;
%             else
%                 hollow = false;
%             end
%         end
    end
    %% overload
    methods
        % we define * call .* as a single mupliy function
        % there * is kron
        function C = mtimes(A,B)
            % only care structure
            if isa(A,'HollowKnight') && isa(B,'HollowKnight')
                if length(A) == 1 && length(B) == 1
                    C = A.*B;
                    return
                end
                % Define mtimes is very difficult. hard coding
                if isrow(A) && isrow(B)
                    C = A(1,1).*B(1,1);
                    for j = 2:size(B,2)
                        D = A(1).*B(j);
                        if ~D.hollow
                            C = [C,D];
                        end
                    end
                    for i = 2:size(A,2)
                        for j = 1:size(B,2)
                            D = A(i).*B(j);
                            if ~D.hollow
                                C = [C,D];
                            end
                        end
                    end
                    C = contractrow(C);
                    return;
                end
                % 1*col = col multipliy means 
                if isrow(A) && ~isrow(B)
                    C = A*B(1,:);
                    for  j= 2:size(B,1)
                        C = [C;A*B(j,:)];
                    end
                    return;
                end
                % col*col  = col*col multipliy means kron
                if ~isrow(A) && ~isrow(B)
                    C = A(1,:)*B;
                    for i = 2:size(A,1)
                        C = [C;A(i,:)*B];
                    end
                    return;
                end
            elseif ~isa(A,'HollowKnight') && isa(B,'HollowKnight')
                % Define mtimes is very difficult. hard coding
                if length(A) == 1 && length(B) == 1
                    C = A.*B;
                    return
                end
                % 1 *(m*n)
                if isscalar(A)
                    C = B;
                    for i = 1:numel(C)
                        C(i) = A.*B(i);
                    end
                    C = contract(C);
                    return
                end
                % row * row
                if isrow(A) && isrow(B)
                    if length(A) == length(B)
                        C = A(1).*B(1);
                        for j = 2:size(B,2)
                            D = A(j).*B(j);
                            if ~D.hollow
                                C = [C,D];
                            end
                        end
                        C = contractrow(C);
                        return;
                    else
                        error('not imply yet')
                    end
                end
                % row*col = 1*? multipliy means
                if isrow(A) && ~isrow(B)
                    if size(A,2) == size(B,1)
                        C = A(1)*B(1,:);
                        for  j= 2:size(B,1)
                            C = [C,A(j)*B(j,:)];
                        end
                        C = contractrow(C);
                    else
                        error('not imply yet')
                    end
                    return;
                end
                % (row*col)*col  = row*? multipliy means kron
                if ~isrow(A) && ~isrow(B)
                    if size(A,2) == size(B,1)
                        C = A(1,:)*B;
                        for i = 2:size(A,1)
                            C = [C;A(i,:)*B];
                        end
                    else
                        error('not imply yet')
                    end
                    return;
                end
            elseif isa(A,'HollowKnight') && ~isa(B,'HollowKnight')
                C = mtimes(B,A);
            else
            end
        end
        function C = times(A,B)
            if isa(A,'HollowKnight') && isa(B,'HollowKnight')
                % return times
                C = innertimes(A,B);
            elseif ~isa(A,'HollowKnight') && isa(B,'HollowKnight')
                C = B;
                if length(A) == 1
                    for i =1:numel(B)
                        C(i).coe = A.*C(i).coe ;
                    end
                elseif isvector(A) && length(A) == size(B,1)
                    for i = 1:size(B,1)
                        C(i,:) = A(i).*B(i,:);
                    end
                elseif isequal(size(A) ,size(B))
                    for i = 1:size(B,1)
                        for j = 1:size(B,2)
                            C(i,j) = A(i,j).*B(i,j);
                        end
                    end
                else
                    error('times wrong')
                end
            elseif isa(A,'HollowKnight') && ~isa(B,'HollowKnight')
                C = times(B,A);% suppose commute!!
            else

            end
        end
        function C = innertimes(A,B)
            C = A               ;
            C.coe = A.coe*B.coe ;
        end
    end
    %% rotation
    methods
        function Am = rotate(A,rotm,rightorleft,options)
            % Input arguments
            %   A         HollowKnight      (bandlimit+1)^2 x 1 complex array, spherical
            %                           harmonics coefficients.
            %   rotm    double      4 x 1 array or 3 x 3 array, rotation in quaternion
            %                           representation or in rotation matrix
            %                           representation, respectively.
            %
            arguments
                A HollowKnight;
                rotm {Spin.mustBeSize(rotm,[3 3;1 5;5 1;1 4;4 1])}= diag([1 1 1]);% [alpha beta gamma]
                rightorleft = 'right';
                options.sym = true;
                options.conjugate = false;
                options.antisymmetry = false;
            end
            optionsCell = namedargs2cell(options);
            Am = rotaterow(A(1,:),rotm,rightorleft,optionsCell{:});
            for i = 2:size(A,1)
                Am = [Am;rotaterow(A(i,:),rotm,rightorleft,optionsCell{:})];
            end
        end
        function A_Lj = rotaterow(A,rotm,rightorleft,options)
            arguments
                A HollowKnight;
                rotm {Spin.mustBeSize(rotm,[3 3;1 5;5 1;1 4;4 1])}= diag([1 1 1]);% [alpha beta gamma]
                rightorleft = 'right';
                options.sym = true;
                options.conjugate = false;
                options.antisymmetry = false;
            end
            optionsCell = namedargs2cell(options);
            A_Lj = contractrow(rotatesingle(A(1),rotm,rightorleft,optionsCell{:}));
            for i = 2:size(A,2)
                Ai = rotatesingle(A(i),rotm,rightorleft,optionsCell{:});
                %                 for k =1:length(Ai)
                %                     if ~isequal(Ai(k).coe,zeros(1,1,class(Ai(k).coe)))
                %                         A_Lj =[A_Lj,Ai(k)];
                %                     end
                %                 end
                A_Lj = [A_Lj,contractrow(Ai)];
            end
            A_Lj = (A_Lj);
        end
        function Ak = rotatesingle(A,rotm,rightorleft,options)
            arguments
                A HollowKnight;
                rotm {Spin.mustBeSize(rotm,[3 3;1 5;5 1;1 4;4 1])}= diag([1 1 1]);% [alpha beta gamma]
                rightorleft = 'right';
                options.sym = true;
                options.conjugate = false;
                options.antisymmetry = false;
            end
            optionsCell = namedargs2cell(options);
            if strcmp(rightorleft,'right')
                RightorLeft = 1; % Check which is right!!!!!
            else
                RightorLeft = -1;
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
                % note the robot toolbos use left axis! clockwise rotation
                % here we change our right rotm to left rotm
                rotm = inv(rotm);
                abc = Oper.Rotation2eul(rotm);% alpha beta gamma in ZYZ
            elseif isequal(size(rotm),[4 1]) || isequal(size(rotm),[1 4])
                if sym(abs(rotm(end))) ~= sym(1)
                    warning('wrong input,eular angle ZYZ right format:[alpha beta gamma det()]');
                    abc = Oper.axang2eul(rotm(1:4));
                    immproper = false;
                else
                    if sym(rotm(end)) == sym(-1)
                        immproper = true;
                    else
                        immproper = false;
                    end
                    abc = rotm(1:3);
                end
            elseif isequal(size(rotm),[5 1]) || isequal(size(rotm),[1 5])
                if sym(abs(rotm(end))) ~= sym(1)
                    error('wrong input,axis angle det right format:[nx ny nz theta det()]');
                end
                if sym(rotm(end)) == sym(-1)
                    immproper = true;
                else
                    immproper = false;
                end
                abc = Oper.axang2eul(rotm(1:4));
            end
            if options.sym
                abc = sym(abc);
            else
                %abc
            end
            Ak = rotateinner(A,abc,RightorLeft,immproper,options.conjugate,options.antisymmetry);
        end
    
    end
    methods(Abstract)
        Ak = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
    end
    %% contract
    methods
        function A = horzcat(varargin)
            arguments(Repeating)
                varargin HollowKnight;
            end
            %             ncolB = size(B,2);
            for i =1:nargin
                classL(i) =  string(class(varargin{i}));
            end
            classname = string(class(varargin{1}));
            chooseL = find(classname ~= classL);
            % force class transformation    
            for i = 1:chooseL
                varargin{i} = feval(classname,varargin{i});
            end
            A = builtin('horzcat',varargin{:});
        end
        function A = vertcat(varargin)
            arguments(Repeating)
                varargin HollowKnight;
            end
            %             ncolB = size(B,2);
            for i =1:nargin
                ncol(i) =  size(varargin{i},2);
            end
            ncolmax = max(ncol);
            %ncolB = size(B,2);
            %ncolC = size(C,2);
            for i = 1:nargin
                if ncolmax ~= ncol(i)
                    % ?
                    varargin{i} = varargin{i}.padded(ncolmax);
                end
            end
            A = builtin('vertcat',varargin{:});
        end
        function A = padded(A,ncol,nrow)
            arguments
                A HollowKnight;
                ncol =1;
                nrow =1;
            end
            if ncol-size(A,2) <= 0
                return;
            end
            nrow = size(A,1);
            A = [A,HollowKnight.creator(ncol-size(A,2),nrow,class(A))];
        end
        function A = cleanrow(A)
            if isrow(A)
                coeList = [A.coe];
                coeList(isnan(coeList)) = 0;
                A(sym(coeList) == sym(0)) = [];
            else
                error('Not a row!')
            end
        end
        function B = contract(A,options)
            arguments
                A HollowKnight;
                options.forgetcoe = false;
            end
            optionsCell = namedargs2cell(options);
            nrow = size(A,1);
            B = contractrow(A(1,:));
            for i = 2:nrow
                %HollowKnightRow = cleanrow();
                HollowKnightRow = contractrow(cleanrow(A(i,:)),optionsCell{:});
                B = [B;HollowKnightRow];
            end
        end
    end
    methods(Abstract)
        B = contractrow(A);
    end
    methods
        function A = HollowMe(A)
            A.coe= nan;
        end
    end
    methods(Static)
        function [comparerow_unique,sumrow_unique] = generalcontractrow(comparerow,sumrow)
            arguments
                comparerow
                sumrow
            end
            % rm nan
            selecL =  ~logical(sum(isnan(sumrow),2));
            %
            [comparerow_unique,sumrow_unique] = HollowKnight.generalcontractrow2(comparerow(selecL,:),sumrow(selecL,:));
            Lic = find(sum(logical(zeros(1,size(sumrow,2),class(sumrow))==sumrow_unique),2));
            sumrow_unique(Lic,:) = [];
            comparerow_unique(Lic,:) = [];
        end
        function [comparerow_unique,sumrow_unique] = generalcontractrow2(comparerow,sumrow)
            arguments
                comparerow
                sumrow
            end
            if size(comparerow,1) == 1
                comparerow_unique = comparerow;
                sumrow_unique = sumrow;
                return;
            end
            [comparerow_unique,~,~] = unique(comparerow,'rows');
            if size(comparerow_unique,1) == size(comparerow,1)
                comparerow_unique = comparerow;
                sumrow_unique = sumrow;
                return;
            end
            sumrow_unique = zeros([size(comparerow_unique,1),size(sumrow,2)],class(sumrow));
            for i = 1:size(comparerow_unique,1)
                comparerow_unique_tmp = comparerow_unique(i,:);
                [~,index_ic] = ismember(comparerow,comparerow_unique_tmp,'rows');
                sumrow_unique(i,:) =  sum(sumrow(logical(index_ic),:),1);
            end
        end
        function HollowKnightObj = creator(ncol,nrow,classname)
            HollowObj =  feval(classname);
            HollowObj = HollowObj.HollowMe();
            HollowKnightObj = repmat(HollowObj,[nrow,ncol]);
        end
        function varargout = StanderdInput(varargin,options)
            arguments(Repeating)
                varargin
            end
            arguments
                options.nrow = 1;
                options.ncol = 1;
            end
            % auto  nQnums nQnums_col
            nrow = options.nrow;
            ncol = options.ncol;
            %
            for i =1:numel(varargin)
                if length(varargin{i}) ==1
                    tmp = repmat(varargin{i},[nrow,ncol]);
                else
                    tmp = varargin{i};
                end
                if iscolumn(tmp)
                    tmp = repmat(tmp,[1,ncol]);
                elseif isrow(tmp)
                    tmp = repmat(tmp,[nrow,1]);
                end
                varargout{i} = tmp;
            end
        end
    end
    %% operator
    methods
        % InnerProduct
        function Smat = InnerProduct(A,B,options)
            arguments
                A HollowKnight
                B HollowKnight
                options.sym = false;
                options.strict = false;
                options.union = false;
                options.square = true;
            end
            optionsCell = namedargs2cell(options);
            if size(A,1)~=size(B,1)
                % try to contract
                A = contract(A);
                B = contract(B);
                if length(A)~=length(B) && options.union
                    % try to union

                end
            end
            % nA = nB>
            if options.sym
                Ssingle = sym(0);
            else
                Ssingle = 0;
            end
            Smat = repmat(Ssingle,[size(A,1),size(B,1)]);
            for i =1:size(A,1)
                for j = 1:size(B,1)
                    Smat(i,j) = InnerProduct_row(A(i,:),B(j,:),optionsCell{:});
                end
            end
        end
        function SingleSum = InnerProduct_row(A_row,B_row,options)
            arguments
                A_row HollowKnight
                B_row HollowKnight
                options.sym = true;
                options.strict = false;
                options.union = false;
                options.square = true;
            end
            if options.sym
                SingleSum = sym(0);
            else
                SingleSum = 0;
            end
            if isrow(A_row) && isrow(B_row)
                for i =1:size(A_row,2)
                    iA_row = A_row(i);
                    if iA_row.hollow
                        continue;
                    end
                    for j = 1:size(B_row,2)
                        jB_row = B_row(j);
                        if jB_row.hollow
                            continue;
                        end
                        if iA_row == jB_row
                            % Wrong:-> We should use Bj for col ; Ai row ?
                            SingleSum = SingleSum + A_row(i).coe*B_row(j).coe;
                        end
                    end
                end
            end
        end
        function HollowKnightObj = Prow(HollowKnightObj,P,options)
            arguments
                HollowKnightObj
                P
                options.antisymmetric = true;
            end
            if isvector(P)
                permVec = P;
                %permMat = HollowKnight.permVecToMat(P);
            else
                permVec = HollowKnight.permMatToVec(P);
                %permMat = P;
            end
            % here we dont consider supersymmetry; if need, the whole class
            % need to be recoded
            if  options.antisymmetric
                Parity = HollowKnight.ParityPerm(permVec)^(HollowKnightObj(1).J*2);
            else
                Parity = 1;
            end
            HollowKnightObj = Parity.*HollowKnightObj(:,permVec);
        end
        function HollowKnightObj = Pcol(HollowKnightObj,P,options)
            arguments
                HollowKnightObj
                P
                options.antisymmetric = false;
            end
            if isvector(P)
                permVec = P;
                %permMat = HollowKnight.permVecToMat(P);
            else
                permVec = HollowKnight.permMatToVec(P);
                %permMat = P;
            end
            % here we dont consider supersymmetry; if need, the whole class
            % need to be recoded
            if  options.antisymmetric
                Parity = HollowKnight.ParityPerm(permVec)^(HollowKnightObj(1).J*2);
            else
                Parity = 1;
            end
            HollowKnightObj = Parity.*HollowKnightObj(permVec,:);
        end
        function HollowKnightObj = Pmat(HollowKnightObj,P,options)
            arguments
                HollowKnightObj HollowKnight;
                P
                options.antisymmetric = false;
            end
            if isvector(P)
                %permVec = P;
                permMat = HollowKnight.permVecToMat(P);
            else
                %permVec = HollowKnight.permMatToVec(P);
                permMat = P;
            end
            HollowKnightObj = permMat * HollowKnightObj;
            %             HollowKnightObj = contractrow(HollowKnightObj(:).');
            %             for i = 2:size(permMat,1)
            %                 HollowKnightObjTmp = permMat(i,:).*HollowKnightObj;
            %                 HollowKnightObjTmp = contractrow(HollowKnightObjTmp(:).');
            %                 HollowKnightObj = [HollowKnightObj;HollowKnightObjTmp];
            %             end
        end
    end
    methods(Static)
        % perm
        function permMat =  permVecToMat(permVec)
            n = length(permVec);
            permMat = zeros(n,n);
            for i = 1:n
                permMat(i, permVec(i)) = 1;
            end
        end
        function permVec =  permMatToVec(permMat)
            n = length(permMat);
            permVec = (1:n)*permMat.';
        end
        function Parity = ParityPerm(permVec)
            if isvector(permVec)
            else
                permVec = Spin.permMatToVec(permVec);
            end
            y = eye(numel(permVec));
            Parity = det( y(:,permVec));
        end
    end
     
end