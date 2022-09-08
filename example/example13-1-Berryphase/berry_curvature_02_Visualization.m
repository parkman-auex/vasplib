%% berry_curvature_02_Visualization
syms h theta phi real;
%
% % phi   = 0;
% sigma_x = [0 1;1 0];
% sigma_y = [0 -1i;1i 0];
% sigma_z = [1 0;0 -1];
% %
% h_vector(1) = h*sin(theta)*cos(phi);
% h_vector(2) = h*sin(theta)*sin(phi);
% h_vector(3) = h*cos(theta);
% 
% H = sigma_x*h_vector(1)+sigma_y*h_vector(2)+sigma_z*h_vector(3);
%%
BC = 0.5*sin(theta);
BC_fun = matlabFunction(BC,'Vars',[h,theta,phi]);
%%
ntheta = 30;
nphi = 60;
nh = 3;
S = zeros(ntheta*nphi*nh,3);
A_S = S;
count = 0;
for itheta = linspace(0,pi,ntheta)
    for phi = linspace(0,2*pi,nphi)
        for h = linspace(1,3,nh)
            count = count+1;
            S(count,:) = [h,itheta,phi];
            A_h  = BC_fun(h,itheta,phi);
            A_S(count,:) = [A_h itheta phi];
        end
    end
end
R =  vasplib.Sph2Cart(S);
A_R = vasplib.Sph2Cart(A_S);

import vasplib_tool.*
vasplib_tool.quiverplot(R,A_R);

