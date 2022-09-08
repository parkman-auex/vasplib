function H00_H11_cell_list = H00_H11_cell_list_gen(H_hr,klist_s,Rm,fin_dir,principle_layer,mode)
if nargin < 6
    mode = 't';
end

switch mode
    case 't'
        disp('hr_dat Hxyz wannierTOOls(wt) ');
        if nargin <5
            principle_layer = 3;
        end        
        if nargin <4
            fin_dir = 3;
        end
        if nargin <3
            POSCAR_read;
        end
        if nargin <2
            [~,~,klist_s,~,~]=kpathgen3D(Rm);
        end
        WAN_NUM = length(H_hr(1).Hnum);
        k_n=length(klist_s(:,1));
        H00_H11_cell_list{1,k_n} = zeros(WAN_NUM); % H00
        H00_H11_cell_list{2,k_n} = zeros(WAN_NUM); % H01
        hr_sparse = hr2hr_sparse(H_hr);
        for j = 1:k_n
            kpoints_f = klist_s(j,:);
            [H_hrz,~,~]= hrz_gen(hr_sparse,kpoints_f,fin_dir);
            [H00,H01,~] = Green_prepare(H_hrz,principle_layer,fin_dir);
            H00_H11_cell_list{1,j} = H00;
            H00_H11_cell_list{2,j} = H01;
        end
    case 'ts'
        
end
     
end

