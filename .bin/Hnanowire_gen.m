function [H_wire,Hnum_list_wire,vector_list_wire] = Hnanowire_gen(Hnum_list,vector_list,Nslab,np)
    %--------  reshape  -------- 
    fin_dir = [3,2,1];
    [vector_list_init,sort_label] = sortrows(vector_list,fin_dir) ;% sort fin_dir
    WAN_NUM=length(Hnum_list{1});
    NRPTS=length(Hnum_list);
    Hnum_list_init = reshape(full(cell2mat(Hnum_list)),WAN_NUM,WAN_NUM,NRPTS);
    Hnum_list_init = Hnum_list_init(:,:,sort_label); 
    % --------  init  --------
    WAN_NUM_x = WAN_NUM*Nslab(1);
    WAN_NUM_y = WAN_NUM_x*Nslab(2);
    % --------  3D begining  --------
    [unique_z,unique_label_z]= unique(vector_list_init(:,3),'rows');
    cutlist_z = unique_label2cutlist(unique_label_z,NRPTS);
    NRPTS_z = length(unique_z);
    vector_list_wire = zeros(NRPTS_z,3);

    % --------  vector  --------
    for iz = 1:NRPTS_z
        vector_list_wire(iz,3) = unique_z(iz);
        vertor_list_xy{iz} = vector_list_init(cutlist_z(iz,1):cutlist_z(iz,2),:);
        Hnum_list_xy{iz}   = Hnum_list_init(:,:,cutlist_z(iz,1):cutlist_z(iz,2));
        Hnum_list_wire{iz} = sparse(WAN_NUM_y,WAN_NUM_y);
    end

    if np >1    
        parfor iz = 1:NRPTS_z
            fprintf('Gen (%d/%d) NRPT z \n',iz,NRPTS_z);
            Hnum_list_xy_iz = Hnum_list_xy{iz};
            vertor_list_xy_iz = vertor_list_xy{iz};
            Hnum_list_wire{iz} = Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WAN_NUM,WAN_NUM_x,WAN_NUM_y,Nslab);     
        end
    else
        for iz = 1:NRPTS_z
            Hnum_list_xy_iz = Hnum_list_xy{iz};
            vertor_list_xy_iz = vertor_list_xy{iz};
            Hnum_list_wire{iz} = Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WAN_NUM,WAN_NUM_x,WAN_NUM_y,Nslab);
        end
    end
    H_wire.HnumL = Hnum_list_wire;
    H_wire.vectorL = vector_list_wire;
end
function Hnum_list_wire_iz =  Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WAN_NUM,WAN_NUM_x,WAN_NUM_y,Nslab)
Hnum_list_wire_iz = sparse(WAN_NUM_y,WAN_NUM_y);
[NRPTS_xy,~] = size(vertor_list_xy_iz);
% --------  2D begining  --------
[unique_y,unique_label_y]= unique(vertor_list_xy_iz(:,2),'rows');
cutlist_y = unique_label2cutlist(unique_label_y,NRPTS_xy);
NRPTS_y = length(unique_y);
% vertor_list_y = repmat([0 0 unique_z(iz)],3,1);
%Hnum_list_iy_iz{NRPTS_y} = sparses(WAN_NUM_x,WAN_NUM_x);
    for iy = 1:NRPTS_y
        fprintf('%d th NRPT z ---- Gen (%d/%d) NRPT y \n',iz,iy,NRPTS_y);
%         vertor_list_y(iy,2) = unique_y(iy);
        vertor_list_x_iy_iz = vertor_list_xy_iz(cutlist_y(iy,1):cutlist_y(iy,2),:);
        Hnum_list_x_iy_iz   = Hnum_list_xy_iz(:,:,cutlist_y(iy,1):cutlist_y(iy,2));
        [NRPTS_x,~] = size(vertor_list_x_iy_iz);
        
        % --------  1D begining  --------
        Hnum_list_iy_iz =sparse(WAN_NUM_x,WAN_NUM_x);
        for ix = 1:NRPTS_x % go over all NRPTS
            fprintf('%d th NRPT z ---- %d th NRPT y ---- Gen (%d/%d) NRPT x \n',iz,iy,ix,NRPTS_x);
            % lattice vector of the hopping
            ind_R_x = vertor_list_x_iy_iz(ix,:);
            jump_fin_x=ind_R_x(1);
            % store by how many cells is the hopping in finite direction
            %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
            % speed up more and more !!
            for  i = 1:WAN_NUM
                for j = 1:WAN_NUM
                    % amplitude of the hop is the same
                    amp = Hnum_list_x_iy_iz(i,j,ix);
                    if norm(amp) > 0
                        for icur_sc_vec = 1:Nslab(1) % go over all cells in finite direction mini
                            hi= i + (icur_sc_vec-1)*WAN_NUM ;
                            %disp(hi);
                            hj= j + (icur_sc_vec+jump_fin_x-1)*WAN_NUM ;
                            %disp(hj);
                            % decide whether this hopping should be added or not
                            to_add=1;
                            %disp(hj);
                            % if edges are not glued then neglect all jumps that spill out
                            if hj <= 0 || hj > WAN_NUM_x
                                to_add=0;
                            end
                            if to_add == 1
                                %Hnum_list_iy_iz(hi,hj) = Hnum_list_iy_iz(hi,hj)+ amp;
                                Hnum_list_iy_iz(hi,hj) = amp;
                            end
                        end
                    end
                end
            end
            %OUT_Hnum_list(:,ih) = temp_Hnum;
        end
        jump_fin_y =  unique_y(iy);
        %           ind_R_y = vertor_list_x_iy_iz(iy,:);        
        for  i = 1:WAN_NUM_x
            for j = 1:WAN_NUM_x
                amp = Hnum_list_iy_iz(i,j);
                if norm(amp) > 0
                    for icur_sc_vec = 1:Nslab(2) % go over all cells in finite direction mini
                        hi= i + (icur_sc_vec-1)*WAN_NUM_x ;
                        %disp(hi);
                        hj= j + (icur_sc_vec+jump_fin_y-1)*WAN_NUM_x ;
                        %disp(hj);
                        % decide whether this hopping should be added or not
                        to_add=1;
                        %disp(hj);
                        % if edges are not glued then neglect all jumps that spill out
                        if hj <= 0 || hj > WAN_NUM_y
                            to_add=0;
                        end
                        if to_add == 1
                            %Hnum_list_wire_iz(hi,hj) = Hnum_list_wire_iz(hi,hj)+ amp;
                            Hnum_list_wire_iz(hi,hj) =  amp;
                        end
                    end
                end
            end
        end
    end
end
function cutlist= unique_label2cutlist(unique_label,NRPTS)
    cutlist(:,1)= unique_label;
    cutlist(1:end-1,2)= unique_label(2:end)-1;
    cutlist(end,2) = NRPTS;
end