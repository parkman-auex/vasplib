%% Cd3As2_type2_formation_low_kp_model
lattice_mode = 'T';
syms A M D F k real;
syms B B_ B3 B3_ k_plus k_minus;
syms E0k real;
syms C0 C1 C2 real;
syms M0 M1 M2 real;
syms B1 B2 kzc real;
syms K_x K_y K_z K_x2 K_y2 K_z2 real;

H_updn=[0        ,A*k_plus,0        ,        0;...
        A*k_minus,0       ,0     ,         0;...
        0        ,0       ,0        ,-A*k_minus;...
        0        ,0       ,-A*k_plus,0         ];

% H_updn=[0        ,A*k_plus,B3_        ,        0;...
%         A*k_minus,0       ,0       ,         B_;...
%         B        , 0      ,0        ,-A*k_minus;...
%           0      ,B3      ,-A*k_plus,0         ];
H_diag=diag([M,-M,M,-M]);
H_E0k=diag([E0k,E0k,E0k,E0k]);
H_gamma=H_E0k+H_updn+H_diag;

disp(H_gamma);


%% note: k_alpha must write into K_x formation

% B =B1  ;
% B_=B1  ;
%  B       = B1*kzc*(K_x2-K_y2+2*1i*K_x*K_y)+B2*kzc*(K_x2-K_y2-2*1i*K_x*K_y);
% B_      = B1*kzc*(K_x2-K_y2-2*1i*K_x*K_y)+B2*kzc*(K_x2-K_y2+2*1i*K_x*K_y);
%   B       = B1*(K_x^2-K_y^2);
%   B_      = B1'*(K_x^2-K_y^2);
%   B3       = B1*(K_x^2-K_y^2);
%   B3_      = B1'*(K_x^2-K_y^2);
%   
% 

%       Ba      =  1*B1*(K_x*K_y);  Bb      =  1i*B2*(K_x*K_y);
%       Bd      =  1*B2*(K_x*K_y);  Bc      =   -1i*B1*(K_x*K_y);


%      Bd     =  1*B2*(K_x^2);  Bc      =  0;
%      Ba      =  0; Bb      =  -1*B2*(K_y^2);

     Ba      =  1*B1*(K_x^2+K_y^2);  Bb      =  B2*(K_x^2-K_y^2);
     Bd      = B2*(K_x^2-K_y^2); Bc      =  -1*B1*(K_x^2+K_y^2);
    
%     Ba      =  1i*B1*(k_plus);  Bb      =  -1i*B2*k_minus;
%     Bd      =  1i* B2*(k_plus); Bc      =  -1* B1*k_plus;
    
%    Ba      =  1i*B1   ;     Bb      =  -1i*B2;
%    Bd      =  1* B2; Bc      =  -1* B1;
% ¡¾B1 B2;B2 B1¡¿
Bk = [Ba, Bb;...
      Bd, Bc];


E0k     = C1*K_z^2;

k_plus  = K_x + 1i* K_y;
k_minus = K_x - 1i* K_y;
M       = M0  -M2*(K_x^2+K_y^2)-M1*(K_z^2);


B_SOTI  = [ [0 0 ;0 0],Bk';Bk,[0 0 ;0 0]];
H_gamma = H_gamma + B_SOTI ;
H=subs(H_gamma);
disp(H);
% eig(H) solve
