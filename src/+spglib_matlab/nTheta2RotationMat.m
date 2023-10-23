function RotationMat = nTheta2RotationMat(n,theta)
n_x = n(1);
n_y = n(2);
n_z = n(3);
%% left
% RotationMat(1,1)=n_x^2*(1-cosd(theta))+cosd(theta);
% RotationMat(1,2)=n_x*n_y*(1-cosd(theta))+n_z*sind(theta);
% RotationMat(1,3)=n_x*n_z*(1-cosd(theta))-n_y*sind(theta);
% RotationMat(2,1)=n_y*n_z*(1-cosd(theta))-n_z*sind(theta);
% RotationMat(2,2)=n_y^2*(1-cosd(theta))+cosd(theta);
% RotationMat(2,3)=n_y*n_z*(1-cosd(theta))+n_x*sind(theta);
% RotationMat(3,1)=n_z*n_x*(1-cosd(theta))+n_y*sind(theta);
% RotationMat(3,2)=n_z*n_y*(1-cosd(theta))-n_x*sind(theta);
% RotationMat(3,3)=n_z^2*(1-cosd(theta))+cosd(theta);
%% right
RotationMat(1,1)=n_x^2*(1-cosd(theta))+cosd(theta);
RotationMat(1,2)=n_x*n_y*(1-cosd(theta))-n_z*sind(theta);
RotationMat(1,3)=n_x*n_z*(1-cosd(theta))+n_y*sind(theta);
RotationMat(2,1)=n_y*n_z*(1-cosd(theta))+n_z*sind(theta);
RotationMat(2,2)=n_y^2*(1-cosd(theta))+cosd(theta);
RotationMat(2,3)=n_y*n_z*(1-cosd(theta))-n_x*sind(theta);
RotationMat(3,1)=n_z*n_x*(1-cosd(theta))-n_y*sind(theta);
RotationMat(3,2)=n_z*n_y*(1-cosd(theta))+n_x*sind(theta);
RotationMat(3,3)=n_z^2*(1-cosd(theta))+cosd(theta);
end