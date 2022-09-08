%% direct kp model
syms k_x k_y k_z a b c real
 K_plus  = k_x + 1i* k_y;
 K_minus = k_x - 1i* k_y;
 K_plus_2 =(k_x + 1i* k_y)^2;
 K_minus_2 =(k_x - 1i* k_y)^2;
 K_x2y2 = k_x^2 + k_y^2;
 K_pm_3 = 2*k_x*(k_x^2 - 3*k_y^2);
 K_x=k_x;
 K_y=k_y;
 K_z=k_z;
 K_x2=k_x^2;
 K_y2=k_y^2;
 K_z2=k_z^2;
 K_1 = k_x;
 K_2 =(k_x + sqrt(3)*k_y)/2;
 K_12 =(k_x - sqrt(3)*k_y)/2;