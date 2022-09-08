%% autohermi_hr  
function H_xyz = autohermi_hr(H_xyz,mode)
if nargin<2
    mode = 'sym';
end
 H_xyz_tmp = H_xyz;
if strcmp(mode,'sym')
    disp('For sym Hr, The requirement is strick! The sym vasiable is real or in the original Hr,use conj()');
   [NRPTS,~]=size(H_xyz);
   disp(NRPTS);
    WAN_NUM=length(H_xyz(1).Hcoe);
    WAN_NUM_half = WAN_NUM;
    V = [H_xyz.vector];
    V =reshape(V,3,length(V)/3)';
    for i = 1:NRPTS
        
        vector_tmp = V(i,:);
        vector_tmp_oppo = -vector_tmp;
        [~,j]=ismember(vector_tmp_oppo,V,'rows');
        fprintf('check %d th : NRPT, the vector is %d %d %d, the opposite vector is %d,%d,%d the %d th NRPT\n',i,...
                                                vector_tmp(1),vector_tmp(2),vector_tmp(3),...
                                                 vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3),j) ;
        if i == j
            if ~isequal(H_xyz(i).Hcoe ,H_xyz_tmp(j).Hcoe')
                fprintf('The homecell hamilton is not hermi, never mind, we will hermi it enforcely!\n');
                disp(H_xyz(i).Hcoe);
                H_xyz_tmp = set_hop(H_xyz(i).Hcoe'/2+H_xyz(j).Hcoe/2,1,1,[0,0,0],H_xyz_tmp,'symmat');
                disp('change into');
                disp(H_xyz_tmp(i).Hcoe);
            end
            continue;
        end
        
        if j == 0
            fprintf('The opposite vector hamilton is not exist, build it!\n');
            H_xyz_tmp = set_hop(H_xyz(i).Hcoe',1,1,-vector_tmp,H_xyz_tmp,'symmat');
            disp(H_xyz_tmp(i).Hcoe');
            V = [V;-vector_tmp];
            continue;
        elseif ~isequal(H_xyz_tmp(i).Hcoe ,H_xyz_tmp(j).Hcoe')
            fprintf('The opposite vector hamilton is not hermi, replace it by strong mat! \n');
            N1 = nnz(H_xyz_tmp(i).Hcoe);N2 = nnz(H_xyz_tmp(j).Hcoe');
            %disp([N1,N2]);
            if N1 >= N2
                fprintf('The %d th NRPT is stonger!\n',i);
                disp(H_xyz_tmp(i).Hcoe);
                H_xyz_tmp = set_hop(H_xyz_tmp(i).Hcoe',1,1,-vector_tmp,H_xyz_tmp,'symmat');
            elseif N1 < N2
                fprintf('The %d th NRPT is stonger!\n',j);
                disp(H_xyz_tmp(j).Hcoe');
                H_xyz_tmp = set_hop(H_xyz_tmp(j).Hcoe',1,1,vector_tmp,H_xyz_tmp,'symmat');
            else
            end
        elseif isequal(H_xyz_tmp(i).Hcoe ,H_xyz_tmp(j).Hcoe')
%             disp(H_xyz(i).Hcoe);
%             disp(H_xyz(j).Hcoe);
            disp('hermi test pasts');
        else
            warning('!!!!!');
        end
        
        
    end
    
    
    
elseif  strcmp(mode,'num')
    
end
  

H_xyz= H_xyz_tmp;
end