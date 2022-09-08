function Parity_list = Parity_eight(Occupi,list)
if nargin <2
    list = [1:8];
end
    for i = 1:length(list)
        Parity_list(i) = Parity('label_file_'+string(list(i)),Occupi);
%         Parity_list.M1  = Parity('label_file_2',Occupi);
%         Parity_list.M2  = Parity('label_file_3',Occupi);
%         Parity_list.M3  = Parity('label_file_4',Occupi);
%         Parity_list.A  = Parity('label_file_5',Occupi);
%         Parity_list.H1  = Parity('label_file_6',Occupi);
%         Parity_list.H2  = Parity('label_file_7',Occupi);
%         Parity_list.H3  = Parity('label_file_8',Occupi);
    end
end
