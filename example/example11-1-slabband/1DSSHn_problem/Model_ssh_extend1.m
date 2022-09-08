%Model_ssh

syms v w v2 w2 u1 u2 real

H00 = [0 v;conj(v) 0];
H01 = [u1 v2;w u2];
H02 = [0 0;w2 0];

SSH_base = init_hr(2);

% set hop
SSH_base = set_hop(H00,0,0,[0 0 0],SSH_base,'symmat');
SSH_base = set_hop(H01,0,0,[1 0 0],SSH_base,'symmat');
SSH_base = set_hop(H02,0,0,[2 0 0],SSH_base,'symmat');
% auto hermi
SSH_extend = autohermi_hr(SSH_base );

% 
% 
% syms a b c k_x k_y k_z real;
% 
% Rm = [a 0 0;0 b 0;0 0 c];
% SSH_reveal = reveal_hr(SSH_extend,Rm);
% 
% 
% % Green function
% u1 = 0;
% u2 = 0;
% a  = 1 ;
% Rm_green = subs(Rm);
% SSH_Green_H = Subsall(SSH_extend,'sym'); 
% SSH_reveal_Green_H = reveal_hr(SSH_Green_H,Rm_green);
% 
% syms w ;
% SSH_Green =inv(w*eye(2)-SSH_reveal_Green_H);