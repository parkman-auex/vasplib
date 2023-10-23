function [Hue,surf_level,hing_level] = orbone2hsv(orb_one,discrimination,center,orientation)
%--------  nargin  --------
if nargin < 2
    discrimination = 0.1;
end
if nargin <3
    center  = [0.5, 0.5,0.5];
end

if nargin <4
    orientation = 3;
end

orb_init = orb_one - center;
%--------  init  --------
switch orientation
    case 1
        x = orb_init(2);
        y = orb_init(3); 
    case 2
        x = orb_init(3);
        y = orb_init(1);
    case 3
        x = orb_init(1);
        y = orb_init(2);
end
r = norm([x y]);
z = x+1i*y;
Hue = (angle(z)+pi)/(2*pi);
G_hinge = ((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
%G_hinge = (r+1i*discrimination-0.5*1.414)^-1;
G_surf = (r+1i*discrimination-0.5)^-1;
hing_level = -imag(G_hinge);
surf_level = -imag(G_surf);
    
end

% function creatHSV_cir(size)
% 
% % Hue Begin
% x_center = ceil(size/2);
% y_center = x_center;
% for m = 1:size
%     for n = 1:size
%         if m<y_center && n>=x_center
%             theta = asin((n-x_center)/sqrt((m-y_center)^2+(n-x_center)^2));
%             else if m>=y_center && n>x_center
%                 theta = asin((m-y_center)/sqrt((m-y_center)^2+(n-x_center)^2))+pi/2;
%                 else if m>y_center && n<=x_center
%                     theta = asin((x_center-n)/sqrt((m-y_center)^2+(n-x_center)^2))+pi;
%                     else if m<=y_center && n<x_center
%                         theta = asin((y_center-m)/sqrt((m-y_center)^2+(n-x_center)^2))+3*pi/2;
%                     end
%                 end
%             end
%         end
%         H(m,n) = theta/(2*pi);
%     end
% end
% H(x_center,y_center) = 0;
% 
% % Hue End
% 
% % Saturation Begin
% x_center = ceil(size/2);
% y_center = x_center;
% for m = 1:size
%     for n = 1:size
%         S(m,n) = (sqrt((m-y_center)^2+(n-x_center)^2))/(size/2);
%     end
% end
% 
% % Saturation End
% V = ones(size);
% HSV = cat(3,H,S,V);
% RGB = hsv2rgb(HSV);
% RGB = cut(RGB,(size/2));
% imshow(RGB);
% 
% function T = cut(T,r)
% Tsize = size(T,1);
% center = ceil(Tsize/2);
% for m =1:Tsize
%     for n =1:Tsize
%         if sqrt((m-center)^2+(n-center)^2)>=r
%             T(m,n,:) = 1;
%         end
%     end
% end
