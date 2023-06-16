function H_spin = add_nn(H_spin, exchange_num, exchange_label)
% add normal exchange terms, like 1nn 2nn in a 2D lattice
arguments
    H_spin SpinModel
    exchange_num int8 = 2
    exchange_label {mustBeA(exchange_label,{'sym','string'})} = ["J1", "J2"]
end
%% check 1
if length(exchange_label) ~= exchange_num
    error("The length of exchange_label must be equal to exchange_num");
end
hr = H_spin.HR_obj_base;
%% check 2
exchange_label = sym(exchange_label,'real');
allvars = [exchange_label, H_spin.exchange_label];

if H_spin.exchange_num + exchange_num ~= length(unique(allvars))
    error("The labels of exchange terms must be individual");
end
%% search bonds, here more likely the intra-layer bonds
search_range = [exchange_num-1, exchange_num-1, 0];  % a proper guess
Accuracy = 1e-2;
Rlength_cut = 12; % a large value for metal system
hr = hr.nn(search_range, Accuracy, Rlength_cut);
%% construct the pseudo-TBmodel
hr = hr.H_TBSK_gen('level_cut', exchange_num, 'SymHopping_Human', exchange_label);
%% check and output to screen
nn_found = length(hr.symvar_list);
if nn_found ~= exchange_num
    disp(hr.nn_information());
    disp(hr.symvar_list);
    error("The number of exchanges found is not equal to input, "+...
        "please check the nn_information listed above");               
end
%% save properties
H_spin.HR_obj = H_spin.HR_obj + hr;
end