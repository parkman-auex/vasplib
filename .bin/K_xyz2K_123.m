%K_xyz2K_123
SQRT3 = sqrt(3);
iSQRT3 = 1i*sqrt(3);
if lattice_mode == 'T'
% cubic
syms K1 K2 K12 K3 K1_2 K2_2 K12_2 K3_2 real
K_x2=K1_2;
K_y2=K2_2;
K_z2=K3_2;
K_x=K1;
K_y=K2;
K_z=K3;
H=subs(H);
disp(H);
end

if lattice_mode == 'T_TB'
% cubic
syms K1 K2 K12 K3 K1_2 K2_2 K12_2 K3_2 real

coskx=1-K1_2/2;
cosky=1-K2_2/2;
coskz=1-K3_2/2;
sinkx=K1;
sinky=K2;
sinkz=K3;
H=subs(H);
disp(H);
end

% % Hexagonal lattice
% if lattice_mode == 'H'
% % K_plus K_minus
% syms K1 K2 K12 K3 K1_2 K2_2 K12_2 K3_2 real
% K_x2y2 = 2/3 * ( K1_2 + K12_2 + K2_2 ) ;
% 
% K_plus = 2/3 * ( K1 + 1/2 * (1+iSQRT3)*K12 + 1/2 * (1 - iSQRT3) * K2 ) ;
% K_minus = 2/3 * ( K1 + 1/2 * (1-iSQRT3)*K12 + 1/2 * (1 + iSQRT3) * K2 );
% 
% K_plus_2 = 2/3 * ( 2 * K1_2 + (-1 - iSQRT3) * K12_2 + ( -1 + iSQRT3 ) * K2_2 ) ;
% K_minus_2 = 2/3 * ( 2 * K1_2 + (-1 + iSQRT3) * K12_2 + ( -1 - iSQRT3 ) * K2_2 );
% 
% K_pm_3 = -K1 + K2 + K12;
% K_z = K3;
% K_z2 = K3_2;
% 
% H=subs(H);
% disp(H);
% end

% Hexagonal lattice 60 
if lattice_mode == 'H'
% K_plus K_minus
syms K1 K2 K12 K3 K1_2 K2_2 K12_2 K3_2 real
K_x2y2 = 2/3 * ( K1_2 + K12_2 + K2_2 ) ;
K_x2my2 = 2/3 * ( 2*K1_2 - K12_2 - K2_2 ) ;

K_plus = 2/3 * ( K1 + exp(pi*1i/3)*K12 + exp(-1i*pi/3) * K2 ) ;
K_minus = 2/3 * ( K1 + exp(-1i*pi/3)*K12 + exp(pi*1i/3) * K2 );

K_plus_2 = 4/3 * ( K1_2 + exp(-2*pi*1i/3)* K12_2 + exp(2*pi*1i/3)  * K2_2);
K_minus_2 = 4/3 * (K1_2 + exp(2*1i*pi/3) * K12_2 + exp(-2*pi*1i/3) * K2_2) ;

K_pm_3 = -K1 + K2 + K12;
K_z = K3;
K_z2 = K3_2;

H=subs(H);
disp(H);
end


% % Hexagonal lattice 120
% if lattice_mode == 'H2'
% % K_plus K_minus
% syms K1 K2 K12 K1p2 K3 K1_2 K2_2 K12_2 K1p2_2 K3_2  real
% K_x2y2 =  2/3 * ( K1_2 + K12_2 + K2_2 ) ;
% K_x2my2 = 2/3 * ( -K1_2 - K12_2 + 2*K2_2 ) ;
% 
% K_plus = 2/3  *( exp(-pi*1i/3) * K1 + exp(pi*1i/3) * K12 +  K2 ) ;
% K_minus = 2/3  *( exp(pi*1i/3) * K1 + exp(-pi*1i/3) * K12 +  K2 ) ;
% 
% K_plus_2  = 4/3  *( exp(2*pi*1i/3) * K1_2 + exp(-2*pi*1i/3) * K12_2  +  K2_2 ) ;
% K_minus_2 = 4/3 *( exp(-2*pi*1i/3) * K1_2 + exp(2*pi*1i/3) * K12_2 +  K2_2 ) ;
% 
% K_pm_3 = -K1 + K2 + K12;
% K_z = K3;
% K_z2 = K3_2;
% 
% H=subs(H);
% disp(H);
% end


% Hexagonal lattice 120
if lattice_mode == 'H2'
% K_plus K_minus
syms K1 K2 K12 K1p2 K3 K1_2 K2_2 K12_2 K1p2_2 K3_2  real
K_x2y2 = 2/3 * ( K1_2 + K1p2_2 + K2_2 ) ;

K_x2my2 = 2/3 * ( K1_2 - K1p2_2 - 2*K2_2 ) ;

K_plus  = 2/3  *( exp(-pi*1i/3) * K1 + exp(pi*1i/3) * K1p2 +  K2 ) ;
K_minus = 2/3  *( exp(pi*1i/3) * K1 + exp(-pi*1i/3) * K1p2 +  K2 ) ;

K_plus_2  = 4/3 * ( exp( 2*pi*1i/3)  * K1_2  + exp(-2*pi*1i/3) * K1p2_2  +  K2_2 ) ;
K_minus_2 = 4/3 * ( exp(-2*pi*1i/3)  * K1_2  + exp( 2*pi*1i/3) * K1p2_2  +  K2_2 ) ;

K_pm_3 = -K1 + K2 + K12;
K_z = K3;
K_z2 = K3_2;

H=subs(H);
disp(H);
end

% syms K1 K2 K3 K1_2 K2_2 K3_2 real
% K_x = (2*K1+K2)/sqrt(3);
% K_y = K2;
% K_z = K3;
% K_x2= K_x^2;
% K_y2= K_y^2;
% K_z2= K_z^2;
% 
% H=subs(H);
% disp(H);