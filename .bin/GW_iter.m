function [Gwl,Gwb,Gwr] = GW_iter(H00,H01,w,eta,mu_max,infinity_small)
%% 
if nargin <6
   infinity_small= 1e-10;
end
%%
if nargin <5
   mu_max = 100;
end
epsilon0 = H00;
alpha0 = H01;
beta0 = H01';

%oumega_bk = E_w;
oumega =  (w+1i*eta)*eye(length(H00));
%save('Green_prepare.mat');
epsilon_s0 = epsilon0;
epsilon_s_bar0 = epsilon0;
%[epsilon,epsilon_s,epsilon_s_bar,~,~] = GW_iter_coff(mu,oumega,alpha0,beta0,epsilon0);

for i = 1:mu_max
    Ml = (oumega - epsilon0)\alpha0;
    Mr = (oumega - epsilon0)\beta0;
    M1 = alpha0 * (Ml);
    M2 = beta0  * (Mr);
    M3 = alpha0 * (Mr);
    M4 = beta0  * (Ml);

    difference = M3+ M4;
    if norm(difference,'fro') < infinity_small
        %disp('reach_accuracy')
        epsilon0 = epsilon0 +M3+ M4;
        epsilon_s0 = epsilon_s0 + M3;
        epsilon_s_bar0 = epsilon_s_bar0 + M4;
        break;
    end
%     disp(max( difference(:)));
    alpha0 = M1;
    beta0 = M2;
    epsilon0 = epsilon0 +M3+ M4;
    epsilon_s0 = epsilon_s0 + M3;
    epsilon_s_bar0 = epsilon_s_bar0 + M4;
    
%     [epsilon,epsilon_s,epsilon_s_bar,alpha,beta,difference] = GW_iter_coff_once(oumega,epsilon0,epsilon_s0,epsilon_s_bar0,alpha0,beta0);
%     alpha0 = alpha;
%     beta0 = beta;
%     epsilon0 = epsilon;
%     epsilon_s0 = epsilon_s;
%     epsilon_s_bar0 = epsilon_s_bar;
    %disp(i);
    %disp(max(difference(:)));

end
% epsilon = epsilon0;
% epsilon_s = epsilon_s0 ;
% epsilon_s_bar = epsilon_s_bar0;
Gwl = inv((oumega-epsilon_s0));
Gwb = inv((oumega-epsilon0));
Gwr = inv((oumega-epsilon_s_bar0));
end
% 
% function [epsilon,epsilon_s,epsilon_s_bar,alpha,beta] = GW_iter_coff(mu,oumega,alpha0,beta0,epsilon0)
% % if nargin <2
% %     disp('final');
% % end
% global infinity_small;
% 
% % iter finish one
% if alpha0 < infinity_small & beta0 < infinity_small
%     disp('temp_return');
%     return;
% end
% 
% % disp([mu max(abs(alpha0(:))) max(abs(beta0(:)))]);
% 
% if mu ==0
%     %load('Green_prepare.mat');
%     alpha = alpha0;
%     beta = beta0;
%     epsilon = epsilon0;
%     epsilon_s = epsilon0;
%     epsilon_s_bar = epsilon0;
%     %mu = mu+1;
% % elseif mu ==1
% %     [epsilon0,epsilon_s0,epsilon_s_bar0,alpha0,beta0,~] = GW_iter_coff(0);  
% % elseif mu>1
% %     [epsilon0,epsilon_s0,epsilon_s_bar0,alpha0,beta0,~] = GW_iter_coff(mu-1,epsilon,epsilon_s,alpha,beta);   
% % end
% elseif mu >=  1
%     [epsilon_,epsilon_s_,epsilon_s_bar_,alpha_,beta_] = GW_iter_coff(mu-1,oumega,alpha0,beta0,epsilon0);
%     [epsilon,epsilon_s,epsilon_s_bar,alpha,beta,difference] = GW_iter_coff_once(oumega,epsilon_,epsilon_s_,epsilon_s_bar_,alpha_,beta_);
%     disp(mu);
%     disp(max(difference(:)));
% end
% end


function [epsilon,epsilon_s,epsilon_s_bar,alpha,beta,difference] = GW_iter_coff_once(oumega,epsilon0,epsilon_s0,epsilon_s_bar0,alpha0,beta0)
%     M1 = alpha0 * inv(oumega - epsilon0)*alpha0;
%     M2 = beta0  * inv(oumega - epsilon0)*beta0;
%     M3 = alpha0 * inv(oumega - epsilon0)*beta0;
%     M4 = beta0  * inv(oumega - epsilon0)*alpha0;
    M1 = alpha0 * ((oumega - epsilon0)\alpha0);
    M2 = beta0  * ((oumega - epsilon0)\beta0);
    M3 = alpha0 * ((oumega - epsilon0)\beta0);
    M4 = beta0  * ((oumega - epsilon0)\alpha0);
    alpha = M1;
    beta = M2;
    epsilon = epsilon0 +M3+ M4;
    difference = epsilon-epsilon0;
%     disp(max( difference(:)));
    epsilon_s = epsilon_s0 + M3;
    epsilon_s_bar = epsilon_s_bar0 + M4;
end