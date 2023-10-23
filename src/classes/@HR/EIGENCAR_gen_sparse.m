function EIGENCARout = EIGENCAR_gen_sparse(H_hr,fermi,norb_enforce,klist_s_tmp)
            % -------------- nargin ------------------
            if nargin <2
                fermi = 0;
            end
            if nargin <3
                norb_enforce  = -1;
            end
            if nargin <4
                if isempty(H_hr.klist_s)
                    H_hr = H_hr.kpathgen3D('KPOINTS');
                end
                klist_s_tmp = H_hr.klist_s;
            end
            % -------------- nargin ------------------
            %disp("EIGENCAR gen for H_xyz(wt TB) type: HR class ");
            vectorlist = double(H_hr.vectorL) ;
            WANNUM = H_hr.WAN_NUM;
            if WANNUM > 500
                print_mode = 1;
            else
                print_mode = 0;
            end
            NRPTS_tmp = H_hr.NRPTS;
            [kn,~] = size(klist_s_tmp);
            %--------  check  --------
            if norb_enforce <0
                NBANDS=WANNUM;
            elseif norb_enforce >0
                NBANDS=norb_enforce;
            else

            end
            EIGENCAR = zeros(NBANDS,kn);
            if kn >1
                for ki =1:kn
                    kx=klist_s_tmp(ki,1);
                    ky=klist_s_tmp(ki,2);
                    kz=klist_s_tmp(ki,3);

                    if strcmp(H_hr.Type,'sparse')
                        Htemp=sparse(WANNUM ,WANNUM);
                        for i=1:NRPTS_tmp
                            Htemp = Htemp +H_hr.HnumL{i}*exp(1i*2*pi*([kx ky kz]*vectorlist(i,:)'));
                        end
                        Hout = Htemp;
                    end
                    if norb_enforce <0
                        [~, U]=eig(full(Hout));
                    elseif norb_enforce >0
                        [~, U]=eigs(Hout,NBANDS,fermi);
                        [~, U]= HR.sorteig(U);
                    else
                    end
                    EIGENCAR(:,ki) = diag(U);
                    if print_mode ==1
                        fprintf('%d th kpoints has been calculated in %d kpoints total\n',ki,kn);
                    end
                end
            else
                kx=klist_s_tmp(1);
                ky=klist_s_tmp(2);
                kz=klist_s_tmp(3);

                if strcmp(H_hr.Type,'sparse')
                    Htemp=sparse(WANNUM ,WANNUM);
                    for i=1:NRPTS_tmp
                        Htemp = Htemp +H_hr.HnumL{i}*exp(1i*2*pi*([kx ky kz]*vectorlist(i,:)'));
                    end
                    Hout = Htemp;
                end
                if norb_enforce <0
                    [A, U]=eig(Hout);
                elseif norb_enforce >0
                    [A, U]=eigs(Hout,NBANDS,fermi);
                    [A, U]= HR.sorteig(U,A);
                else
                end
                EIGENCAR = diag(U);
                WAVECAR = A;
                if print_mode ==1
                    fprintf('only %d th kpoint(%7.5f %7.5f %7.5f) has been calculated in %d kpoints total\n',1,kx,ky,kz,kn);
                end
            end
            if kn >1
                EIGENCARout = EIGENCAR;
                %EIGENCARout{2} = WAVECAR;
            else
                EIGENCARout.EIGENCAR= EIGENCAR;
                EIGENCARout.WAVECAR = WAVECAR;
            end
        end
        