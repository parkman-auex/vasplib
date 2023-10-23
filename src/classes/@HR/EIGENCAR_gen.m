function varargout = EIGENCAR_gen(H_hr,options)
arguments
    H_hr HR;
    options.fermi double = 0;
    options.norb double = -1;
    options.klist double = H_hr.klist_frac;
    options.para  = [];
    options.paraname ;
    options.convention {mustBeMember(options.convention,{'I','II'})}= 'II';
    options.show = false;
    options.ax = handle([]);
    options.WEIGHTCAR = false;
    options.ProjectionMethod {mustBeMember(options.ProjectionMethod,{'hinge','surf','select','select-points','slab'})}= 'hinge';
    options.ProjectionStruct = struct('discrimination',0.1,'center',repmat(0.5,[1 H_hr.Dim]),'orientation',3,'sign',true);
    options.printmode = true;
    options.LWAVE = true;
    options.Umat = [];
    options.Oper = [];
    options.subband = [];
    options.returnH = [];
    options.Hermite = true;
end
% -------------- define ------------------
norb_enforce  = options.norb;
Hermite = options.Hermite;
NRPTS_tmp = H_hr.NRPTS;
% -------------- plot ------------------
if options.show
    if  isempty(options.ax)
        ax = BZplot(vasplibobj.Rm,'color','r');
    else
        ax = options.ax;
    end
end
% -------------- nargin ------------------
if length(options.fermi) == 1
    fermi = options.fermi;
    dynamic_fermi = false;
else
    fermi_list = options.fermi;
    dynamic_fermi = true;
end
if isempty(options.klist)
    H_hr = H_hr.kpathgen3D('KPOINTS');
    klist_frac_tmp = H_hr.klist_frac;
    klist_cart = klist_frac_tmp*H_hr.Gk;
else
    klist_frac_tmp = options.klist;
    klist_cart = klist_frac_tmp*H_hr.Gk;
end
if isempty(options.Umat)
    UmatMode = false;
else
    UmatMode = true;
    Umat = options.Umat;
    if size(Umat,3) >1
        LDynamicU = true;
    else
        LDynamicU = false;
    end
    if isempty(options.subband)
        subband = 1:H_hr.WAN_NUM/2;
    else
        subband = options.subband;
    end
end
if options.WEIGHTCAR
    orb_list  = H_hr.orbL;
    signlist =1;
    %HSVCAR_hinge = vasplib.HSVCAR_gen(orb_list,'hinge',0.05,[0.5,0.5,0.5],-3);
    HSVCAR = vasplib.HSVCAR_gen(orb_list,options.ProjectionMethod,...
        options.ProjectionStruct.discrimination,...
        options.ProjectionStruct.center,...
        options.ProjectionStruct.orientation...
        );
    switch options.ProjectionMethod
        case 'hinge'
            
        case 'surf'
      
        case {'select','select-points'}
            signlist = sign(HSVCAR(:,1));
        case 'slab'
            if options.ProjectionStruct.sign
                switch  options.ProjectionStruct.orientation
                    case 1
                        signlist = sign((orb_list(:,1)-0.5));
                    case 2
                        signlist = sign((orb_list(:,2)-0.5));
                    case 3
                        signlist = sign((orb_list(:,3)-0.5));
                    case 4
                        signlist = sign((orb_list(:,4)-0.5));
                end
                signlist(HSVCAR(:,1) == 0) = 0;
            end
    end
end
if isempty(options.Oper)
    OperMode = 0;
else
    % options.WEIGHTCAR = true;
    OperMode = 1;
end
% -------------- nargin ------------------
if H_hr.overlap
    NRPTS_tmp_S = size(H_hr.vectorL_overlap,1);
end
if strcmp(H_hr.Type,'list')
    if H_hr.overlap
        H_hr = H_hr.SliceGen();
        [ij_list_S,index_row_S] = sortrows(H_hr.vectorL_overlap(:,4:5));
        H_hr = H_hr.reseq(':',index_row,index_row_S);
        [Sparse_vector_S,SliceList_S] = unique(ij_list_S,'rows');
        N_Sparse_vector_S = size(Sparse_vector_S,1);
        CutList_S = [SliceList_S,[SliceList_S(2:end)-1;NRPTS_tmp_S]];
    else
        H_hr = H_hr.SliceGen();
    end
end
%disp("EIGENCAR gen for H_xyz(wt TB) type: HR class ");
% reseq makes HR different! take attention
HnumList = H_hr.HnumL ;
vectorList = double(H_hr.vectorL(:,1:H_hr.Dim)) ;
if H_hr.overlap
    Snum_list = double(H_hr.SnumL) ;
    vectorlist_overlap = double(H_hr.vectorL_overlap(:,1:H_hr.Dim)) ;
end
if UmatMode
    WANNUM = length(subband);
else
    WANNUM = H_hr.WAN_NUM;
end
if options.printmode
    pb = vasplib_tool_outer.CmdLineProgressBar('BAND calculating ');
end
[kn,~] = size(klist_frac_tmp);
%--------  check  --------
if norb_enforce <0
    NBANDS=WANNUM;
elseif norb_enforce >0
    NBANDS=norb_enforce;
else

end
if options.LWAVE
    if isempty(options.subband)
        WAVECAR  = zeros(WANNUM,NBANDS,kn);
        subband = 1:NBANDS;
    else
        subband = options.subband;
        WAVECAR  = zeros(WANNUM,length(subband),kn);
    end
else
    WAVECAR = [];
end
EIGENCAR = zeros(NBANDS,kn);
if options.WEIGHTCAR || OperMode
    WEIGHTCAR = EIGENCAR;
end
% give back
WANNUM = H_hr.WAN_NUM;
% Factorlist_R $H_{i j}^{\mathbf{k}}=\left\langle\chi_{i}^{\mathbf{k}}|H| \chi_{j}^{\mathbf{k}}\right\rangle=\sum_{\mathbf{R}} e^{i \mathbf{k} \cdot\left(\mathbf{R}+\mathbf{t}_{j}-\mathbf{t}_{i}\right)} H_{i j}(\mathbf{R})$
FactorList = exp(1i*2*pi*vectorList*klist_frac_tmp.');
if H_hr.overlap
    FactorListS = exp(1i*2*pi*vectorlist_overlap*klist_frac_tmp.');
end
% convention I
if strcmp(options.convention,'I') && (isempty(H_hr.tjmti) || size(H_hr.tjmti{1},1)~= WANNUM)
    H_hr = H_hr.tjmti_gen();
end
if strcmp(options.convention ,'I')
    tji_mat_r = H_hr.tjmti{1}; % tj - ti
    % efactor tji 
    for i = 1:H_hr.Dim
        kjiL{i} =  pagemtimes(tji_mat_r(:,:,i),reshape(klist_cart(:,i),[1 1 kn]));
    end
    FactorList_tji = exp(1i*(fold(@plus,kjiL)));
end
if strcmp(options.convention,'II')
    for ki =1:kn
        if dynamic_fermi
            fermi = fermi_list(ki);
        end
        FactorListki = FactorList(:,ki);
        %Hmat_tij = FactorList_tij(:,:,ki);
        if strcmp(H_hr.Type,'sparse')
            Htemp=sparse(WANNUM ,WANNUM);
            for i=1:NRPTS_tmp
                Htemp = Htemp + (HnumList{i}*FactorListki(i));
            end
            %Hout = Htemp;
            if Hermite
                Hout = (Htemp+Htemp')/2;
            else
                Hout = Htemp;
            end
            %Hout = sparse(Hout);
        elseif strcmp(H_hr.Type,'list')
            Hout = zeros(WANNUM ,WANNUM);
            Hnum_list_k = HnumList.*FactorListki;
            for i=1:H_hr.N_Sparse_vector
                Hout(H_hr.Sparse_vector(i,1),H_hr.Sparse_vector(i,2)) = sum(Hnum_list_k(H_hr.CutList(i,1):H_hr.CutList(i,2)));
            end
            if Hermite
                Hout = (Hout+Hout')/2;
            end
            %
            if H_hr.overlap
                Snum_list_k = Snum_list.*FactorListS(:,ki);
                Sout = zeros(WANNUM , WANNUM);
                for i=1:N_Sparse_vector_S
                    Sout(Sparse_vector_S(i,1),Sparse_vector_S(i,2)) = sum(Snum_list_k(CutList_S(i,1):CutList_S(i,2)));
                end
                Sout = (Sout+Sout')/2;
                %Hout = (Sout)*Hout*inv(Sout);
            end
        else
            Hout = sum(pagemtimes(HnumList,reshape(FactorListki,[1 1 NRPTS_tmp])),3);
            if Hermite
                Hout = (Hout+Hout')/2;
            end
            if H_hr.overlap
                Sout = sum(pagemtimes(SnumList,reshape(FactorListS(:,ki),[1 1 NRPTS_tmp_S])),3);
                Sout = (Sout+Sout')/2;
            end
        end
        if UmatMode
            Hout = Umat\Hout*Umat;
            Hout = Hout(subband,subband);
            if Hermite
                Hout = (Hout+Hout')/2;
            end
        end
        if norb_enforce <0
            if H_hr.overlap
                [A, U]=eig(full(Hout),full(Sout));
            else
                [A, U]=eig(full(Hout));
            end
        elseif norb_enforce >0
            if H_hr.overlap
                [A, U]=eigs(Hout,Sout,NBANDS,fermi);
                [A, U]= park.sorteig(U,A);
            else
                [A, U]=eigs(Hout,NBANDS,fermi);
                [A, U]= park.sorteig(U,A);
            end
        else
        end
        EIGENCAR(:,ki) = diag(U);
        if options.LWAVE
            WAVECAR(:,:,ki) = A(:,subband);
        end
        if options.WEIGHTCAR 
            [~,WEIGHTCAR(:,ki)] = vasplib.COLORCAR_gen(A,HSVCAR,signlist);
        end
        if OperMode
            [WEIGHTCAR(:,ki)] = vasplib.Observecar_gen(A,options.Oper);
        end
        if options.printmode
            pb.print(ki,kn,' ...');
        end
    end
    % normalize phases to get u instead of phi
    %                 for j =1:size(WAVECAR,3)
    %                     WAVECAR(:,:,j) = WAVECAR(:,:,j).* exp(-2i*pi*(H_hr.orbL*klist_frac_tmp(j,:).'));
    %                 end
elseif strcmp(options.convention,'I')
    tji_mat = double(H_hr.tjmti{2});
    for ki =1:kn
        if dynamic_fermi
            fermi = fermi_list(ki);
        end
        FactorListki = FactorList(:,ki);
        Hmat_tji = FactorList_tji(:,:,ki);
        if strcmp(H_hr.Type,'sparse')
            %tij_phase = exp(1i*2*pi*(tij_mat(:,:,1)*kx+tij_mat(:,:,2)*ky+tij_mat(:,:,3)*kz));
            % debug try
            %tij_mat(:,:,1) = tij_mat(:,:,1).';
            %tij_mat(:,:,2) = tij_mat(:,:,2).';
            %tij_mat(:,:,3) = tij_mat(:,:,3).';
            k1=klist_frac_tmp(ki,1);
            k2=klist_frac_tmp(ki,2);
            k3=klist_frac_tmp(ki,3);
            Htemp=zeros(WANNUM ,WANNUM);
            for i=1:NRPTS_tmp
                Htemp = Htemp +full(HnumList{i})...
                    .* exp(...
                    1i*2*pi*(...
                    (tji_mat(:,:,1)+vectorlist(i,1))*k1 + ...
                    (tji_mat(:,:,2)+vectorlist(i,2))*k2 + ...
                    (tji_mat(:,:,3)+vectorlist(i,3))*k3  ...
                    )...
                    );
            end
            if Hermite
                Hout = (Htemp+Htemp')/2;
            else
                Hout = Htemp;
            end
            %Hout = sparse(Hout);
        elseif strcmp(H_hr.Type,'list')
            Hout = zeros(WANNUM ,WANNUM);
            Hnum_list_k = HnumList.*FactorListki;
            for i=1:H_hr.N_Sparse_vector
                Hout(H_hr.Sparse_vector(i,1),H_hr.Sparse_vector(i,2)) = sum(Hnum_list_k(H_hr.CutList(i,1):H_hr.CutList(i,2)));
            end
            Hout = Hout.*Hmat_tji;
            if Hermite
                Hout = (Hout+Hout')/2;
            end
            %
            if H_hr.overlap
                Snum_list_k = Snum_list.*FactorListS(:,ki);
                Sout = zeros(WANNUM , WANNUM);
                for i=1:N_Sparse_vector_S
                    Sout(Sparse_vector_S(i,1),Sparse_vector_S(i,2)) = sum(Snum_list_k(CutList_S(i,1):CutList_S(i,2)));
                end
                Sout = Sout.*Hmat_tji;
                Sout = (Sout+Sout')/2;
            end
        else
            Hout = sum(pagemtimes(HnumList,reshape(FactorListki,[1 1 NRPTS_tmp])),3).*Hmat_tji;
            if Hermite
                Hout = (Hout+Hout')/2;
            end
            if H_hr.overlap
                Sout = sum(pagemtimes(SnumList,reshape(FactorListS(:,ki),[1 1 NRPTS_tmp_S])),3).*Hmat_tji;% ?
                Sout = (Sout+Sout')/2;
            end
        end
        if UmatMode
            if LDynamicU
                Hout = Umat(:,:,ki)\Hout*Umat(:,:,ki);
            else
                Hout = Umat\Hout*Umat;
            end
            % debug
            % max(max(abs(Hout((length(subband)+1):end,subband))))
            Hout = Hout(subband,subband);
            if Hermite
                Hout = (Hout+Hout')/2;
            end
        end
        if norb_enforce <0
            if H_hr.overlap
                [A, U]=eig(full(Hout),full(Sout));
            else
                [A, U]=eig(full(Hout));
            end
        elseif norb_enforce >0
            if H_hr.overlap
                [A, U]=eigs(Hout,Sout,NBANDS,fermi);
            else
                [A, U]=eigs(Hout,NBANDS,fermi);
                [A, U]= HR.sorteig(U,A);
            end
        else
        end
        EIGENCAR(:,ki) = diag(U);
        if options.LWAVE
            WAVECAR(:,:,ki) = A;
        end
        if options.WEIGHTCAR
            [~,WEIGHTCAR(:,ki)] = vasplib.COLORCAR_gen(A,HSVCAR,signlist);
        end
        if OperMode
            [WEIGHTCAR(:,ki)] = vasplib.Observecar_gen(A,options.Oper);
        end
        if options.printmode
            pb.print(ki,kn,' ...');
        end
    end
else
end
varargout{1} = EIGENCAR ;
varargout{2} = WAVECAR;
if options.WEIGHTCAR  || OperMode
    varargout{3} = WEIGHTCAR;
end
if options.show
%     [varargout{3},varargout{4}] = H_hr.klist_show(...
%         'klist',klist_frac_tmp*H_hr.Gk,...
%         'fig',fig,...
%         'ax',ax);
end
if options.returnH
    varargout{3} = Hout;% for debug
end
if options.printmode
    pb.delete;
end
%             if kn >1
%                 EIGENCARout = EIGENCAR;
%                 %EIGENCARout{2} = WAVECAR;
%             else
%                 EIGENCARout.EIGENCAR= EIGENCAR;
%                 EIGENCARout.WAVECAR = WAVECAR;
%             end
end
