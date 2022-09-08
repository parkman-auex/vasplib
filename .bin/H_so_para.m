%% SOC Hamilton  band  


%% H_s_so

Hs_so=[0 0;
       0 0];


%% H_p_so

%       px  py  pz  px  py   pz
% Hp_so=[  0  1i   0   0   0    1 ;...%px
%        -1i   0   0   0   0   1i ;...%py 
%          0   0   0  -1 -1i    0 ;...%pz 
%          0   0  -1   0  1i    0 ;...%px 
%          0   0  1i  1i   0    0 ;...%py 
%          1 -1i   0   0   0    0 ;...%pz
%       ];
  
  
% Wangzhijun PHYSICAL REVIEW B 88, 125427 (2013) 
%       px  py  pz  px  py   pz
Hp_so=[  0 -1i   0   0   0    1 ;...%px
        1i   0   0   0   0  -1i ;...%py 
         0   0   0  -1  1i    0 ;...%pz 
         0   0  -1   0  1i    0 ;...%px 
         0   0 -1i  -1i  0    0 ;...%py 
         1  1i   0   0   0    0 ;...%pz
      ];



%% H_d_so


