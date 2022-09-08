%% Extrat Information

%% init -------------------------------------------------------------------
%
first_load = 0;
spilt_data = 1;
neckname = "NaTmO2";
SG_166_flat     = 1;
% define U_struct
U_list = [0,1,3,5,7,9]; % which U list
%       -------------------------------------------------------------------

if first_load == 1
    [~,DIR_list_str]= system("find . -name 'OUTCAR' |awk -F '/OUTCAR' '{print $1}' ");
    system("find . -name 'OUTCAR' |awk -F '/OUTCAR' '{print $1}' >DIR_list.dat ");
    DIR_list = strsplit(strtrim(DIR_list_str),'\n');
    % HT_MAG_TABLE = table();
    HT_MAG_TABLE = AutoJ(DIR_list{1});
    for i = 2:length(DIR_list)
        fprintf('>======= %s ========<\n',DIR_list{i});
        TEMP_tabel = AutoJ(DIR_list{i});
        HT_MAG_TABLE = [HT_MAG_TABLE;TEMP_tabel ];
        fprintf('*******************************************\n');
    end
    writetable(HT_MAG_TABLE,'HT_MAG_TABLE.xlsx');
    writetable(HT_MAG_TABLE,'HT_MAG_TABLE.dat');
end


HT_MAG_TABLE = readtable('HT_MAG_TABLE.xlsx');

if spilt_data == 1
    mkdir('data_collection');
    count = 0;
    for U = [0,1,3,5,7,9]
        rows = HT_MAG_TABLE.U_eV == U;
        count = count +1;
        HT_MAG_TABLE_Ulist{count} = HT_MAG_TABLE (rows,:) ;
        writetable(HT_MAG_TABLE_Ulist{count},"data_collection/HT_MAG_TABLE_U_"+string(U)+'.dat');
        for type = ["nosoc_col","nosoc_ncol","soc_ncol"]
            switch type
                case "nosoc_col"
                    rows = HT_MAG_TABLE_Ulist{count}.LNONCOLLINEAR == 0;
                    HT_MAG_TABLE_NOSOC_COL{count}  = HT_MAG_TABLE_Ulist{count}(rows,:) ;
                    writetable(HT_MAG_TABLE_NOSOC_COL{count},"data_collection/HT_MAG_TABLE_U_"+string(U)+"_nosoc_col"+'.dat');
                case "nosoc_ncol"
                    rows = (HT_MAG_TABLE_Ulist{count}.LNONCOLLINEAR == 1 & HT_MAG_TABLE_Ulist{count}.LSORBIT == 0);
                    HT_MAG_TABLE_NOSOC_NCOL{count} = HT_MAG_TABLE_Ulist{count}(rows,:) ;
                    writetable(HT_MAG_TABLE_NOSOC_NCOL{count},"data_collection/HT_MAG_TABLE_U_"+string(U)+"_nosoc_ncol"+'.dat');
                    for dir = ["z","x"]
                        switch dir
                            case 'z'
                                rows = (contains(HT_MAG_TABLE_NOSOC_NCOL{count}.SYSTEM,"z")...
                                    | contains(HT_MAG_TABLE_NOSOC_NCOL{count}.SYSTEM,"up"));
                                HT_MAG_TABLE_NOSOC_NCOL_z{count} = HT_MAG_TABLE_NOSOC_NCOL{count}(rows,:) ;
                                writetable(HT_MAG_TABLE_NOSOC_NCOL_z{count},"data_collection/HT_MAG_TABLE_U_"+string(U)+"_nosoc_ncol_z"+'.dat');
                            case 'x'
                                rows = (~contains(HT_MAG_TABLE_NOSOC_NCOL{count}.SYSTEM,"z") ...
                                    & ~contains(HT_MAG_TABLE_NOSOC_NCOL{count}.SYSTEM,"up"));
                                HT_MAG_TABLE_NOSOC_NCOL_x{count} = HT_MAG_TABLE_NOSOC_NCOL{count}(rows,:) ;
                                writetable(HT_MAG_TABLE_NOSOC_NCOL_x{count},"data_collection/HT_MAG_TABLE_U_"+string(U)+"_nosoc_ncol_x"+'.dat');
                        end
                    end
                case "soc_ncol"
                    rows = HT_MAG_TABLE_Ulist{count}.LSORBIT == 1;
                    HT_MAG_TABLE_SOC_NCOL{count} = HT_MAG_TABLE_Ulist{count}(rows,:) ;
                    writetable(HT_MAG_TABLE_SOC_NCOL{count},"data_collection/HT_MAG_TABLE_U_"+string(U)+"_soc_ncol"+'.dat');
                    for dir = ["z","x"]
                        switch dir
                            case 'z'
                                rows = (contains(HT_MAG_TABLE_SOC_NCOL{count}.SYSTEM,"z")...
                                    | contains(HT_MAG_TABLE_SOC_NCOL{count}.SYSTEM,"up"));
                                HT_MAG_TABLE_SOC_NCOL_z{count} = HT_MAG_TABLE_SOC_NCOL{count}(rows,:) ;
                                writetable(HT_MAG_TABLE_SOC_NCOL_z{count},"data_collection/HT_MAG_TABLE_U_"+string(U)+"_soc_ncol_z"+'.dat');
                            case 'x'
                                rows = (~contains(HT_MAG_TABLE_SOC_NCOL{count}.SYSTEM,"z") ...
                                    & ~contains(HT_MAG_TABLE_SOC_NCOL{count}.SYSTEM,"up"));
                                HT_MAG_TABLE_SOC_NCOL_x{count} = HT_MAG_TABLE_SOC_NCOL{count}(rows,:) ;
                                writetable(HT_MAG_TABLE_SOC_NCOL_x{count},"data_collection/HT_MAG_TABLE_U_"+string(U)+"_soc_ncol_x"+'.dat');
                        end
                    end
            end
        end
    end
end
