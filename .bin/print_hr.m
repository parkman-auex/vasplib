%% Tool function
function Hout = print_hr(H_hr,print_list,mode)
fprintf("usage: print_hr(H_hr,[1,3,4],'Hnum-list')\n");
fprintf("usage: print_hr(H_hr,[1,3,4],'Hnum-sum')\n");
fprintf("Oprional: Hnum-sum,Hnum-list,Hcoe-list,Hcoe-sum\n");
if nargin < 3
    mode = 'Hnum-list';
end
if strcmp(mode,'Hnum-sum')
    Hout = H_hr(print_list(1)).Hnum;
    for i = 2:length(print_list)
        Hout = Hout +H_hr(print_list(i)).Hnum;
    end
elseif strcmp(mode,'Hnum-list')
    %
    for i = 1:length(print_list)
        disp([i,H_hr(print_list(i)).vector]);
        disp(H_hr(print_list(i)).Hnum);
    end
elseif strcmp(mode,'Hcoe-list')
    %
    for i = 1:length(print_list)
        disp([i,H_hr(print_list(i)).vector]);
        disp(H_hr(print_list(i)).Hcoe);
    end
elseif strcmp(mode,'Hcoe-sum')
        Hout = H_hr(print_list(1)).coe;
    for i = 2:length(print_list)
                Hout = Hout +H_hr(print_list(i)).Hcoe;
    end
end
end