function [Hc_cell_list,T_cell_list,H_cell_list,H_list] = Hc_T_split(HcDesignData)
    H = HcDesignData.H;
    Hc = HcDesignData.Hc;
    Hc = Hc;
    T = HcDesignData.T;
    cutlist = find(isnan(H)==1);
    ncutlist = length(cutlist);
    nH = length(H);
    cut_label_list(1,:) = [1,cutlist(1)-1];
    
    for i = 1:ncutlist/2
        if 2*i+1 >ncutlist
            cut2 = nH;
        else
            cut2 =cutlist(2*i+1)-1;
        end
        cut1 = cutlist(2*i)+1;
        cut_label_list(i+1,:) =[cut1,cut2]; 
    end
    [N_Hc_list,~] = size(cut_label_list);
    for i = 1:N_Hc_list
        Hc_cell_list{i} = Hc(cut_label_list(i,1):cut_label_list(i,2));
        T_cell_list{i} = T(cut_label_list(i,1):cut_label_list(i,2));
        H_cell_list{i} = H(cut_label_list(i,1):cut_label_list(i,2));
        H_list(i) = round(mean(H_cell_list{i}));
    end
    
end