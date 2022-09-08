classdef group
    %% An operation of discontinus group.
    %
    % * Label: sym group
    %
    %% Description of the class:
    %     U : double (optional)
    %         The unitary action on the Hilbert space.
    %         May be None, 
    %
    %     * Late 1700s- Joseph-Louis Lagrange (1736-1813) 利用置换的概念，理解了三次和四次方程为什么有解。（Paolo Ruffini利用同样的思想证明了一般的五次或以上方程没有根式解）
    %     * Early 1800s-Évariste Galois (killed in a duel in 1832 at age 20), and Niels
    %     * Abel (died in 1829 at age 26 of TB) 阐述了代数方程的可解性和置换群的联系，真正利用群的概念解决了这个难题。
    %     * 第一个正式的群的概念由Cayley于1854提出。Cayley, Hamilton, and Sylvester 引进了矩阵的表达。1890年Fedorov发展了各种对称性操作的数学方法，证明了晶体有且仅有230种空间群，次年Schönflies也独立的证明了这个结果。1920年群论开始应用于物理、化学等领域。里面还有很多著名数学家和物理学家的贡献，不一一而举。
    properties
        U;
    end
    properties(GetAccess = protected,Hidden = true)
        e = false;
        close = false;
    end
    methods
        function Gi = group(U)
            %构造此类的实例
            %   G = group(U)
            arguments
                U = nan;
            end
            Gi.U = U;
        end
    end
    %% Group Defination
    % * 单位元：∃𝑒∈𝑆,∀𝑥∈𝑆,𝑒⋅𝑥=𝑥⋅𝑒=𝑥;
    % Identity element
    % There exists an element 𝑒⋅𝑥=𝑥⋅𝑒=𝑥; It is called the identity element of the group.
    methods(Static)
        function E = identity()
            E = group();
            E.e = true;
        end
    end
    methods
        % E
        function E = E(group)
            %    Return identity element with the same structure as self
            if isnan(group.U)
            else
                group.U = eye(length(group.U),class(group.U));
            end
            E = group;
            E.e = true;
        end
    end
    %  Associativity
    % For all a, b, c in G, one has (a・b)・c=a・(b・c). Alreay imply in ();
    % Inverse element
    % Close: ∀𝑥,𝑦∈𝑆,𝑥⋅𝑦∈𝑆;  imply in * ;
    %% overload method
    methods
        % ==
        function TrueOrFalse = eq(G1,G2)
            % return first
            %
            m = length(G1);
            n = length(G2);
            if m == 1 && n ==1
                if G1.e && G2.e
                    TrueOrFalse = true;
                    return;
                end
                if isequal(G1.U,G2.U)
                    TrueOrFalse = true;
                else
                    TrueOrFalse = false;
                end
                return;
            elseif m == 0 && n == 0
                TrueOrFalse =true;
                return;
            elseif m == 0 || n == 0
                TrueOrFalse = false;
                return
            elseif m == 1
                TrueOrFalse(n) = false;
                for i = 1:n
                    TrueOrFalse(i) = eq(G1,G2(i));
                end
                return;
            elseif n == 1
                TrueOrFalse(m) = false;
                for i = 1:m
                    TrueOrFalse(i) = eq(G1(i),G2);
                end
                return;
            else
                % unique need sort( sort need overload lt(less than))
                % G1 = unique(G1);
                % G2 = unique(G2);
                %
                % m = length(G1);
                % n = length(G2);
                if m == n
                    TrueOrFalse = false(size(G1));
                    for i =1:m
                        TrueOrFalse(i) = eq(G1(i),G2(i));
                    end
                else
                    % generate_group ?
                    TrueOrFalse =false;
                end
            end
        end
        % ~=
        function TrueOrFalse=ne(G1,G2)
            TrueOrFalse = ~(G1  == G2);
        end
        % <
        function TrueOrFalse = lt(G1,G2)
            % Sort group elements:
            if G1.e && ~G2.e
                TrueOrFalse = true;
                return;
            elseif ~G1.e && G2.e
                TrueOrFalse = false;
                return;
            else
            end
            L1 = real(G1.U(:)) < real(G2.U(:));
            B1 = ~(real(G1.U(:)) == real(G2.U(:)));
            L2 = imag(G1.U(:)) < imag(G2.U(:));
            B2 = ~(imag(G1.U(:)) == imag(G2.U(:)));
            for i =1:length(B1)
                if B1(i)
                    TrueOrFalse = L1(i);
                    return;
                end
            end
            for i =1:length(B2)
                if B2(i)
                    TrueOrFalse = L2(i);
                    return;
                end
            end
            TrueOrFalse  = false;
        end
        % >
        function TrueOrFalse = gt(G1,G2)
            TrueOrFalse = lt(G2,G1);
        end
        % inv
        function G = inv(G1)
            % Invert Operator
            % Standerd definition
            % need reload /
            if length(G1) == 1
                G =G1;
                if isnan(G.U)
                else
                    G.U = inv(G.U);
                end
            else
                error('not support G-1 yet.');
            end
        end
        % uminus
        function G = uminus(G)
            for i = 1:numel(G)
                G(i).U = -G(i).U;
            end
        end
        % +
        function G = plus(G1,G2)
            G = unique([G1,G2]);
        end
        % -
        function G = minus(G1,G2)
            % unique need sort( sort need overload lt(less than))
            G1 = unique(G1);
            G2 = unique(G2);
            for i = 1:numel(G2)
                temp_label = G1 ==G2(i);
                if sum(temp_label)
                    G1(temp_label) = [];
                end
            end
            G = G1;
        end
        % /
        function G = mrdivide(G1,G2)
            % unique need sort( sort need overload lt(less than))
            G1 = unique(G1);
            G2 = unique(G2);
            %
            m = length(G1);
            n = length(G2);
            G = [];
            if mod(m,n) == 0
                G= quotient(G2,G1);
            else
                G= [];
            end
        end
        % ./
        function G3 = rdivide(G1,G2)
            % unique need sort( sort need overload lt(less than))
            G1 = unique(G1);
            G2 = unique(G2);
            %
            m = length(G1);
            n = length(G2);
            if m == 1 && n ==1
                if G1.e
                    G3 = G2;
                    G3.U = inv(G2);
                    return;
                elseif G2.e
                    G3 = G1;
                    return;
                else
                    G3 = G1;
                    G3.U = G1.*G2.inv();
                end
            else
                error('.*/ is elemenary operator');
            end
        end
        % *
        % we define * call .* as a single mupliy function
        function G3 = mtimes(G1,G2)
            % unique need sort( sort need overload lt(less than))
            G1 = unique(G1);
            G2 = unique(G2);
            %
            m = length(G1);
            n = length(G2);
            % need reload .*
            % important
            newgroup= repmat(G1(1),[1 m*n]);
            count =  0;
            % generate new elements by multiplying old elements with generators
            for i = 1:m
                for j = 1:n
                    count =  count + 1;
                    newgroup(count) = G1(i).*G2(j);
                end
            end
            G3 = unique(newgroup);
        end
        % .*
        function G3 = times(G1,G2)
            %
            m = length(G1);
            n = length(G2);
            if m == 1 && n == 1
                if G1.e
                    G3 = G2;
                    G3.U = inv(G2);
                    return;
                elseif G2.e
                    G3 = G1;
                    return;
                else
                    G3 = G1;
                    G3.U = G1.U*G2.U;
                    if isa(G3.U,'sym')
                        G3.U = simplify(G3.U);
                    end
                end
            else
                error('.* is elemenary operator');
            end
        end
        % ^
        function G3 = power(G1,n)
            G3 = G1.E();
            for i = 1:abs(n)
                G3 =G3.*G1;
            end
            if n >=0

            else
                G3 =G3.inv();
            end
        end
        % sort
        function [group_list,indSort] = sort(group_list)
            % Use insertsort 
            [group_list,indSort] = group.insertsort(group_list);
        end
    end
    %% Group theory
    methods
        function order = order(SymOper)
            order = length(SymOper);
        end
    end
    methods
        % Closure
        function TrurOrFalse = closure(group)
            TrurOrFalse = true;
            if all([group.close])
                return
            else
                for i = 1:numel(group)
                    if  find(ismember(group(i)*group,group) == 0)
                        TrurOrFalse = false;
                        return;
                    end
                end
            end
        end
        % Abelian
        function TrurOrFalse = abelian(group)
            TrurOrFalse = true;
            for Gi = group
                for Gj = group
                    if Gi*Gj ~=Gj*Gi
                        TrurOrFalse =false;
                        return;
                    end
                end
            end
        end
        % generator
        function group = generate_group(gens)
            %     Generate group from gens
            %
            %     Parameters
            %     ----------
            %     gens : iterable of  generator of a group
            %
            %     Returns
            %     -------
            %     group 
            %
            % need reload minus
            gens = unique(gens);
            % here we keep all the elements generated so far
            group = gens;
            % these are the elements generated in the previous step
            oldgroup = gens;
            %            fprintf('group muplicity: %d\n',group.order);
            while true
                newgroup= unique(oldgroup * gens);
                %                fprintf('newgroup muplicity: %d\n',newgroup.order);
                % only keep those that are new
                newgroup = newgroup - group;
                %                fprintf('newgroup muplicity (after -): %d\n',newgroup.order);
                % if there are any new, add them to group and set them as old
                if ~isempty(newgroup)
                    group = unique([group,newgroup]);
                    oldgroup = newgroup;
                    % if there were no new elements, we are done
                else
                    break;
                end
            end
            %group = group;
        end
        function generator = generator(group,options)
            arguments
                group
                options.fast = true;
            end
            if ~options.fast
                if ~group.closure
                    group = group.generate_group();
                end
            end
            Order = group.order();
            pool = 1:Order;
            %ngenerator = 1;
            findit = false;
            count = 0;
            for i = pool
                ChooseL = nchoosek(pool,i);
                for j = 1:size(ChooseL,1)
                    jgroup = group(ChooseL(j,:));
                    if isequal(length(jgroup.generate_group),Order)
                        findit = true;
                        count = count+1;
                        generator{count} = jgroup;
                        ngenerator = i;
                    end
                end
                if findit
                    break;
                end
            end
        end
        function Rank = rank(group)
            if ~group.closure
                group = group.generate_group();
            end
            pool = 1:numel(group);
            for i = plist
                ChooseL = nchoosek(pool,i);
                for j = 1:size(ChooseL,1)
                    jgroup = group(ChooseL(j,:));
                    if isequal(jgroup.generate_group,group)
                        Rank = i;
                        break;
                    end
                end
            end
        end
        % subgroup
        function subgroup = subgroup(group)
            %     Generate all subgroups of group, including the trivial group
            %     and itself.
            %
            %     Parameters
            %     ----------
            %     group : set of PointGroupElement
            %         A closed group as set of its elements.
            %
            %     Returns
            %     -------
            %     subgroups : dict of forzenset: set
            %         frozesets are the subgroups, sets are a generator set of the
            %         subgroup.
            %
            %     # Make all the cyclic subgroups generated by one element
            % hard coding , solve it if needed
            count = 0;
            count = count + 1;
            subgroup{1} = group(1).E;
            count = count + 1;
            subgroup{2} = group;
            if isprime(group.order())
                % 𝑝(𝑝是素数)阶群𝐺均是Abel群,且均同构于整数模𝑝的加法群ℤ𝑝.
                %for igroup = (group(subgroup{1} ~= group))
                %    count = count + 1;
                %    subgroup{count} = igroup;
                %end
            else
                % Lagrange定理 G的每个子群的阶数都是𝐺的阶数的因数.
                plist = group.SolvingFactor(group.order());
                plist(1:2) = [];
                pool = 1:numel(group);
                for i = plist
                    ChooseL = nchoosek(pool,i);
                    for j = 1:size(ChooseL,1)
                        igroup = group(ChooseL(j,:));
                        if igroup.closure
                            count = count + 1;
                            subgroup{count} = igroup;
                        end
                    end
                end
            end
        end
        % Normal subgroup
        function NormalSubgroup = normalsubgroup(group)
            count = 0;
            count = count + 1;
            NormalSubgroup{count} = group(1).E;
            count = count + 1;
            NormalSubgroup{count} = group;
            if isprime(group.order())
                % 𝑝(𝑝是素数)阶群𝐺均是Abel群,且均同构于整数模𝑝的加法群ℤ𝑝.
                return;
            end
            % 若G是交换群, 则G的所有子群都是正规子群。
            if group.abelian
                NormalSubgroup = subgroup(group);
            end
            % Lagrange定理 + 正规子群定义
            ConjugationClassifyCollection  = conjugationseparation(group);
            nC = length(ConjugationClassifyCollection);
            pool = 1:nC;
            plist = group.SolvingFactor(group.order());
            plist(1:2) = [];
            for i = 1:nC
                ChooseL = nchoosek(pool,i);
                for j = 1:size(ChooseL,1)
                    ConjugationClassifyCollectionTmp = ConjugationClassifyCollection(ChooseL(j,:));
                    SubgroupTmp = [ConjugationClassifyCollectionTmp{:}];
                    nSi = length(SubgroupTmp);
                    if ismember(nSi,plist)
                        if SubgroupTmp.closure()
                            count = count + 1;
                            NormalSubgroup{count} = SubgroupTmp;
                        end
                    else
                        continue;
                    end
                end
            end
        end
        % coset
        function CosetClassify  = cosetseparation(H,G,leftorright)
            arguments
                H
                G
                leftorright = 'l';
            end
            pool_origin = G-H;
            pool = pool_origin;
            count = 1;
            CosetClassify{1}=H;
            for iG = pool_origin
                if ismember(iG,pool)
                    count =count+1;
                    switch leftorright
                        case 'l'
                            Coset = iG*H;
                        case 'r'
                            Coset = H*iG;
                    end
                    CosetClassify{count}=Coset;
                    pool = pool - Coset;
                end
            end
        end
        % 商群
        function TrueOrFalse = NSG(H,G)
            TrueOrFalse = true;
            if ~H.closure  
                TrueOrFalse = false;
                return;
            end
            pool_origin = H;
            pool = pool_origin;
            for Hi = pool_origin
                if ismember(Hi,pool)
                    ConjugationClass = conjugation(Hi,G);
                    if ~all(ismember(ConjugationClass,H))
                        TrueOrFalse = false;
                        return;
                    end
                    pool = pool - ConjugationClass;
                end
            end
        end
        function factor_group = quotient(H,G,leftorright)
            arguments
                H
                G
                leftorright = 'l';
            end
            if H.NSG(G)
                factor_group = cosetseparation(H,G,leftorright);
            else
                factor_group = [];
            end
        end
        % 共轭
        function TrueOrFalse = conjugate(G1,G2,G)
            TrueOrFalse = false;
            for i = 1:numel(G)
                Gi = G(i);
                if G1 == Gi*G2*Gi.inv()
                    TrueOrFalse = true;
                    return;
                end
            end
        end
        %Conjugation class
        function ConjugationClass  = conjugation(G1,G)
            ConjugationClass = G1;
            for i = 1:numel(G)
                Gi = G(i);
                if conjugate(G1,Gi,G)
                    ConjugationClass = [ConjugationClass,Gi];
                end
            end
            ConjugationClass = unique(ConjugationClass);
        end
        %Conjugation classify
        function ConjugationClassify  = conjugationseparation(G)
            ConjugationClassify{1} = G;
            ConjugationClassify{1}(:) = [];
            count = 0;
            pool_origin = G;
            pool = pool_origin;
            for Gi = pool_origin
                if ismember(Gi,pool)
                    count = count +1 ;
                    ConjugationClassify{count} = conjugation(Gi,G);
                    pool = pool - ConjugationClassify{count};
                end
            end
        end
        % 类
    end
    %% tools
    methods(Static)
        % algorithm
        function [sortA,indSortA] = insertsort(A)
            len = length(A);
            if isrow(A)
                indSortA = 1:len;
            else
                indSortA = 1:len;
                indSortA = indSortA.';
            end
            for w = 1:len-1
                % 从余下序列中取出一个元素插入到新序列中
                for v = w+1:-1:2
                    try
                        tmplogical = logical(A(v)<A(v-1)) ;
                    catch
                        tmpsym = (A(v)<A(v-1)) ;
                        if isa(tmpsym,'sym')
                            [~,seq_list] = sort([lhs(tmpsym),rhs(tmpsym)]);
                            if seq_list(1)<seq_list(2)
                                tmplogical = true;
                            else
                                tmplogical = false;
                            end
                        else
                            'Unknown error!';
                        end
                    end
                    if tmplogical
                        % 如果新元素小于所插入位置的元素，则交换
                        tmp = A(v-1);
                        tmp_inder = indSortA(v-1);
                        A(v-1) = A(v);
                        indSortA(v-1)=indSortA(v);
                        A(v) = tmp;
                        indSortA(v) =  tmp_inder;
                    else
                        % 否则保持位置不变
                        break;
                    end
                end
            end
            sortA  = A;
        end
        %
        function plist = SolvingFactor(n)
            arguments
                n {mustBeInteger};
            end
            p=2;
            plist(1) = 1;
            plist(2) = n;
            count = 2;
            while p*p <= n
                if mod(n,p) == 0
                    count = count+1;
                    plist(count) = p;
                    count = count+1;
                    plist(count) = n/p;
                end
                p=p+1;
            end
        end
    end
end