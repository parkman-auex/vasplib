%% useful_tools
sigma_0 = [1 0;0 1];sigma_x = [0 1;1 0];sigma_y = [0 -1i;1i 0];sigma_z = [1 0;0 -1];
tau_0 = [1 0;0 1];tau_x = [0 1;1 0];tau_y = [0 -1i;1i 0];tau_z = [1 0;0 -1];
%% Model
syms k_x k_y k_z real;
syms Delta m v lambda;
%
H_2D = (-Delta + m*k_x^2+m*k_y^2)*kron(sigma_z,tau_0)+...
    v*k_x*kron(sigma_x,tau_x)+...
    v*k_y*kron(sigma_x,tau_y)+...
    lambda*(k_x^2-k_y^2)*kron(sigma_z,tau_x)+...
    2*lambda*(k_x * k_y)*kron(sigma_z,tau_y);
%\left(\begin{array}{cccc} m\,{k_{x}}^2+m\,{k_{y}}^2-2\,\Delta  & \lambda \,\left({k_{x}}^2-{k_{y}}^2\right)-k_{x}\,k_{y}\,\lambda \,2{}\mathrm{i} & 0 & k_{x}\,v-k_{y}\,v\,1{}\mathrm{i}\\ \lambda \,\left({k_{x}}^2-{k_{y}}^2\right)+k_{x}\,k_{y}\,\lambda \,2{}\mathrm{i} & m\,{k_{x}}^2+m\,{k_{y}}^2-2\,\Delta  & k_{x}\,v+k_{y}\,v\,1{}\mathrm{i} & 0\\ 0 & k_{x}\,v-k_{y}\,v\,1{}\mathrm{i} & -m\,{k_{x}}^2-m\,{k_{y}}^2+2\,\Delta  & -\lambda \,\left({k_{x}}^2-{k_{y}}^2\right)+k_{x}\,k_{y}\,\lambda \,2{}\mathrm{i}\\ k_{x}\,v+k_{y}\,v\,1{}\mathrm{i} & 0 & -\lambda \,\left({k_{x}}^2-{k_{y}}^2\right)-k_{x}\,k_{y}\,\lambda \,2{}\mathrm{i} & -m\,{k_{x}}^2-m\,{k_{y}}^2+2\,\Delta  \end{array}\right)'
%% S
S = kron(sigma_x,sigma_0)

H_2D_check = S*H_2D*inv(S)




