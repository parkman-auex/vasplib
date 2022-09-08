
 %% prework
 H_tempstring=  string(expand(H));
 H_tempstring=strrep(H_tempstring,'K1^2','K1_2');% for transformation
 H_tempstring=strrep(H_tempstring,'K2^2','K2_2');% for transformation
 H_tempstring=strrep(H_tempstring,'K12^2','K12_2');% for transformation
 H_tempstring=strrep(H_tempstring,'K1p2^2','K1p2_2');% for transformation
 H_tempstring=strrep(H_tempstring,'K3^2','K3_2');% for transformation
 H=str2sym(H_tempstring);
 clear H_tempstring;
 
%% Gamma discrete
% syms k_x k_y k_z a b c real
%  % cubic
%  K_x=(exp(1i*k_x*a)-exp(-1i*k_x*a))/(2*a*1i);
%  K_y=(exp(1i*k_y*b)-exp(-1i*k_y*b))/(2*b*1i);
%  K_z=(exp(1i*k_z*c)-exp(-1i*k_z*c))/(2*c*1i);
%  K_x2=(2-exp(1i*k_x*a)-exp(-1i*k_x*a))/a^2;
%  K_y2=(2-exp(1i*k_y*b)-exp(-1i*k_y*b))/b^2;
%  K_z2=(2-exp(1i*k_z*c)-exp(-1i*k_z*c))/c^2;
 
 % general
 syms k_1p k_12p k_1p2p k_2p k_3p k_1m k_12m  k_1p2m k_2m k_3m a b c real
 
 K1    = (exp(k_1p)-exp(k_1m))/(2*a*1i);
 K12   = (exp(k_12p)-exp(k_12m))/(2*a*1i);
 K1p2   = (exp(k_1p2p)-exp(k_1p2m))/(2*a*1i);
 K2    = (exp(k_2p)-exp(k_2m))/(2*b*1i);
 
 K3    = (exp(k_3p)-exp(k_3m))/(2*c*1i);
 
 K1_2  = (2 - exp(k_1p) - exp(k_1m) ) / a^2;
 K12_2 = (2 - exp(k_12p)- exp(k_12m)) / a^2;
 K1p2_2 = (2 - exp(k_1p2p)- exp(k_1p2m)) / a^2;
 K2_2  = (2 - exp(k_2p) - exp(k_2m) ) / b^2;
 
 K3_2  = (2-exp(k_3p)-exp(k_3m))/c^2;
 
 
  
 H=subs(H);