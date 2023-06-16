function check_exchange_value(H_spin)
arguments
    H_spin SpinModel
end
if length(H_spin.exchange_value) ~= H_spin.exchange_num
     error("The length of exchange_value is inconsistent with exchange_num")
end
disp("The exchange value you set is:")
for i = 1:H_spin.exchange_num
    disp("  " + string(H_spin.exchange_label(i)) + " = " + H_spin.exchange_value(i) + " meV");
end      
end