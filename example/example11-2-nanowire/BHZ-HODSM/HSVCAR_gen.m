function HSVCAR = HSVCAR_gen(orb_list,mode,discrimination,center,orientation)

%--------  nargin  --------
if nargin < 3
    discrimination = 0.1;
end
if nargin < 4
    center  = [0.5, 0.5,0.5];
end

if nargin < 5
    orientation = 3;
end
if nargin < 2
    mode = 'hinge';
end
[norb,~] = size(orb_list);
HSVCAR  = zeros(norb,1);
switch mode
    case 'hinge'
        for i = 1:norb
            orb_one = orb_list(i,:);
            [~,~,HSVCAR(i)] = orbone2hsv(orb_one,discrimination,center,orientation);
        end
    case 'surf'
        for i = 1:norb
            orb_one = orb_list(i,:);
            [~,HSVCAR(i),~] = orbone2hsv(orb_one,discrimination,center,orientation);
        end    
    case 'orient'
        for i = 1:norb
            orb_one = orb_list(i,:);
            [HSVCAR(i),~,~] = orbone2hsv(orb_one,discrimination,center,orientation);
        end
end
% 0-1
WAN_NUM = norb;
HSVCAR = (-normalize(HSVCAR,'range')+1)/2;
HSVCAR(:,2:3) = ones(WAN_NUM,2);
end