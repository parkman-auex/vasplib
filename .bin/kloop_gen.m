function kloop = kloop_gen(H_hk,input_mat,mode)
            import spglib_matlab.*;
            if nargin < 3
                mode = 'kline';
            end
            if nargin < 2
                input_mat = ...
                    [ 0, 0,  0 ;... %kpoints_f
                    0, 0,  1 ;... %fin_dir
                    0, 1,101];%start end nodes
            end
            %%
            kpoints_f = input_mat(1,:);
            if strcmp(mode,'kline')
                fin_dir_list = input_mat(2,:);
                fin_dir = fin_dir_list ==1;
                nodes = input_mat(3,3);
                if nodes <2
                    error('nodes >= 2!!!!');
                end
                k_start = input_mat(3,1);k_end = input_mat(3,2);
                dk = (k_end-k_start)/(nodes-1);
                klists =repmat(kpoints_f,[nodes,1]);
                klists(:,fin_dir) = (k_start:dk:k_end).';
                kloop = klists*H_hk.Gk;
            elseif  strcmp(mode,'kplane')
                kpoints_r = kpoints_f *H_hk.Gk;
                n_vector = input_mat(2,:);
                nx =n_vector(1);ny =n_vector(1);nz = n_vector(3);
                nodes = input_mat(3,3);
                dk = input_mat(3,1);
                dtheta = 360/(nodes-1);
                theta_init_d = input_mat(3,2);
                if norm(abs(n_vector)-[1,0,0])<1e-8
                    disp('rotation along x axis.');
                    dkx = 0;
                    dky = 0;
                    dkz = dk;
                    dk_default_init = [dkx,dky,dkz];
                elseif ny~=0 && nz~=0
                    dkx = 0;
                    dky = sqrt(nz^2/(ny^2+nz^2)) * dk;
                    dkz = nz*dky/ny;
                    dk_default_init = [dkx,dky,dkz];
                elseif ny==0 && nz~=0
                    dkx = 0;
                    dky = dk;
                    dkz = 0;
                    dk_default_init = [dkx,dky,dkz];
                elseif ny~=0 && nz==0
                    dkx = 0;
                    dky = 0;
                    dkz = dk;
                    dk_default_init = [dkx,dky,dkz];
                else
                    dk_default_init = [dk,0,0];
                end
                rotation_mat_init = spglib_matlab.nTheta2RotationMat(n_vector,theta_init_d);
                dk_init = ((rotation_mat_init*dk_default_init.').'*H_hk.Gk).';
                kloop = zeros(nodes,3);
                kloop(1,:) = dk_init.';
                for i =2:nodes
                    rotation_mat = spglib_matlab.nTheta2RotationMat(n_vector,theta_init_d+(i-1)*dtheta);
                    kloop(i,:) =(rotation_mat*dk_init).';
                end
                kloop = kloop+kpoints_r;
            else
                
            end
        end