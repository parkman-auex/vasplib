function J = CacJ(TARGET_TABLE,choose_system_list)
    switch length(choose_system_list) 
        case 4
            response_temp_mat =  [10,11,12];
        case 3
            response_temp_mat =  [10,11];
        case 2
            response_temp_mat =  [10];            
    end
        
    TEMP_TABLE = TARGET_TABLE(choose_system_list,[response_temp_mat,4,3]);
    temp_mat = table2array(TEMP_TABLE);
    J = temp_mat(1:length(choose_system_list),1:length(choose_system_list))\temp_mat(:,end);
    for i =1:length(response_temp_mat)
        fprintf('J%d : %f meV \n',i,J(i)*1000);
    end
end