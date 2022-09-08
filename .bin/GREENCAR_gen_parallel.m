function GREENCAR = GREENCAR_gen_parallel(np,w_list,eta,H00_H01_cell_list,mode)
if nargin < 5
    mode = 'Green_iter';
end
fprintf('You have use %d cores \n',np);
switch mode
    case 'Green_iter'
        %H00_H01_cell_list = H00_H01_cell_list;
        H00_H01_mat = cell2mat(H00_H01_cell_list);
        [WAN_NUM_2,kn] = size(H00_H01_mat);
        WAN_NUM = WAN_NUM_2/2;       
        Nsize1 = length(w_list);
        Nsize2 = kn/WAN_NUM;
        GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        H00_list = reshape(H00_H01_mat(1:WAN_NUM,:),WAN_NUM,WAN_NUM,Nsize2);
        H01_list = reshape(H00_H01_mat(WAN_NUM+1:end,:),WAN_NUM,WAN_NUM,Nsize2);
        parfor i = 1:Nsize1
            w =w_list(i);
            for j= 1:Nsize2
                H00 = H00_list(:,:,j);
                H01 = H01_list(:,:,j);
                [GREENCAR3{i,j},GREENCAR1{i,j},GREENCAR2{i,j}] = GW_iter(H00,H01,w,eta);
            end
        end
        GREENCAR.bulk = GREENCAR1;
        GREENCAR.surf_l = GREENCAR2;
        GREENCAR.surf_r = GREENCAR3;
    case 'Tmatrix_iter'
        %H00_H01_cell_list = H00_H01_cell_list;
        H00_H01_mat = cell2mat(H00_H01_cell_list);
        [WAN_NUM_2,kn] = size(H00_H01_mat);
        WAN_NUM = WAN_NUM_2/2;       
        Nsize1 = length(w_list);
        Nsize2 = kn/WAN_NUM;
        GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        H00_list = reshape(H00_H01_mat(1:WAN_NUM,:),WAN_NUM,WAN_NUM,Nsize2);
        H01_list = reshape(H00_H01_mat(WAN_NUM+1:end,:),WAN_NUM,WAN_NUM,Nsize2);
        parfor i = 1:Nsize1
            w =w_list(i);
            for j= 1:Nsize2
                H00 = H00_list(:,:,j);
                H01 = H01_list(:,:,j);
                [GREENCAR3{i,j},GREENCAR1{i,j},GREENCAR2{i,j}] = Tmatrix_iter(H00,H01,w,eta);
                % [GREENCAR1{i,j}, GREENCAR2{i,j},GREENCAR3{i,j}]=   Tmatrix2Green00(Tmatrix,H00,H01,w,eta);
            end
        end
        GREENCAR.bulk = GREENCAR1;
        GREENCAR.surf_l = GREENCAR2;
        GREENCAR.surf_r = GREENCAR3;
    case  't'
        if iscell(H00_H01_cell_list)
            disp('do nothing');
        elseif isstruct(H00_H01_cell_list)
            H_hr = H00_H01_cell_list;
            H00_H01_cell_list =H_cell_list_gen(H_hr);
        end
        Nsize1 = length(w_list);
        Nsize2 = length(H00_H01_cell_list);
        %w_list= linspace(w_range(1),w_range(2),w_number);
        WAN_NUM = length(H00_H01_cell_list{1});
        GREENCAR{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        for i = 1:Nsize1
            for j= 1:Nsize2
                GREENCAR{i,j} =  GREENCAR_one_gen(H00_H01_cell_list{j},w_list(i),eta,'TB');
            end
        end
    case 'Tmatrix_eig'
        %H00_H01_cell_list = H00_H01_cell_list;
        WAN_NUM = length(H00_H01_cell_list{1,1});
        [~,kn] =size( H00_H01_cell_list);
        Nsize1 = length(w_list);
        Nsize2 = kn;
        GREENCAR{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        for i = 1:Nsize1
            for j= 1:Nsize2
                w =w_list(i);
                H00 = H00_H01_cell_list{1,j};
                H01 = H00_H01_cell_list{2,j};
                Tmatrix = Tmatrix_gen(H00,H01,w,0);
                GREENCAR{i,j} = Tmatrix2pband(Tmatrix,H00,H01,w);
            end
        end
    case 'Tmatrix_inv'
        H00_H01_mat = cell2mat(H00_H01_cell_list);
        [WAN_NUM_2,kn] = size(H00_H01_mat);
        WAN_NUM = WAN_NUM_2/2;       
        Nsize1 = length(w_list);
        Nsize2 = kn/WAN_NUM;
        GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        H00_list = reshape(H00_H01_mat(1:WAN_NUM,:),WAN_NUM,WAN_NUM,Nsize2);
        H01_list = reshape(H00_H01_mat(WAN_NUM+1:end,:),WAN_NUM,WAN_NUM,Nsize2);
        parfor i = 1:Nsize1
            w =w_list(i);
            for j= 1:Nsize2
                H00 = H00_list(:,:,j);
                H01 = H01_list(:,:,j);             
                Tmatrix = Tmatrix_gen(H00,H01,w,eta);
                [GREENCAR1{i,j}, GREENCAR2{i,j},GREENCAR3{i,j}]=   Tmatrix2Green00(Tmatrix,H00,H01,w,eta);
            end
        end
        GREENCAR.bulk = GREENCAR1;
        GREENCAR.surf_l = GREENCAR2;
        GREENCAR.surf_r = GREENCAR3;


end

% delete(np_handle);
end

