function [M_cell,H_cell,T_cell] = MHT_cell_from_table_struct(table_struct,seq)
    for i = 1:length(seq)
        M_cell{i} = table_struct(seq(i)).datatable.M;
        H_cell{i} = table_struct(seq(i)).datatable.H;
        T_cell{i} = table_struct(seq(i)).datatable.T;
    end
end