function H_spin = add_interlayer(H_spin, exchange_num, exchange_label, Rlength_cut, directions)
% long-range exchange terms, for example, between different layers
arguments
    H_spin SpinModel
    exchange_num int8 = 1
    exchange_label {mustBeA(exchange_label,{'sym','string'})} = "Js"    
    Rlength_cut double = 20; % set a smaller length can speed up
    directions {mustBeMember(directions,{'x', 'y', 'z'})} = 'z';
    %int8 = 3 % 1 2 3 for x y z
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
%% check 3
switch directions
    case 'x'
        directions = 1;
    case 'y'
        directions = 2;
    case 'z'
        directions = 3;
end
%% search long-range bonds
search_range = [1, 1, 1];  % a small guess
Accuracy = 1e-2;
 % a very large value
hr = hr.nn(search_range, Accuracy, Rlength_cut);
%% cut-off the bonds in x-y plane
ns = hr.nn_store;
ns_cut = ns(ns(:,2+directions) ~= 0, :);
ns_cut(:,10) = ns_cut(:,10) - min(ns_cut(:,10)) + 1;
hr.nn_store = ns_cut;
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