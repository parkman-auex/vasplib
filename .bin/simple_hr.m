%%
function [H_hr,line000] = simple_hr(H_hr,mode)
if nargin <2
    mode = 'Hxyz';
end

switch mode
    case 'Hxyz'
        [NRPTS,~]=size(H_hr);
        WAN_NUM=length(H_hr(1).Hnum);
        ZerosMAT = zeros(WAN_NUM);
        rm_list=[];
        for i = 1:NRPTS
            if isequal(H_hr(i).Hnum,ZerosMAT)
                rm_list=[rm_list,i];
            end
            if H_hr(i).Degen >1
                H_hr(i).Hnum  = H_hr(i).Hnum./H_hr(i).Degen ;
                H_hr(i).Degen = 1;
            end
        end
        H_hr(rm_list,:) = [];
        zerovecter = [0 , 0, 0];
        
        V = [H_hr.vector];
        V =reshape(V,3,length(V)/3)';
        % disp(vector);
        % disp(V);
        [~,line000]=ismember(zerovecter,V,'rows');
    case 'sparse'
        disp(" Simple hr for H_xyz(wt TB) type: sparse ");
        Hnum_list = H_hr.HnumL ;
        vectorlist = H_hr.vectorL ;
        [NRPTS,~]=size(vectorlist);
        WAN_NUM = length(Hnum_list{1});
        ZerosMAT = sparse(WAN_NUM,WAN_NUM);
        rm_list=[];
        for i = 1:NRPTS
            if isequal(Hnum_list{i},ZerosMAT)
                rm_list=[rm_list,i];
            end
        end
        Hnum_list{rm_list} = [];
        vectorlist(rm_list,:) =[];
        zerovecter = [0 , 0, 0];
        [~,line000]=ismember(zerovecter,vectorlist,'rows');
        H_hr.HnumL = Hnum_list ;
        H_hr.vectorL = vectorlist;
end


end