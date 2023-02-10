%kmesh_gen generate an uniform kmesh list
% [klist_s, klist_r] = kmesh_gen(opts)
% klist_s : Fractional klist
% klist_r : Cartesian klist (Ang^-1 unit)
%
% see http://www.wanniertools.com/input.html#kcube-bulk for illustration
% mode = corner(default) or center
%   corner mode is same as wanniertools; center mode will move the start
%   point of vk to let the original point in the center of kcube
% edge = half(default) or full
%   "half" means only half points on the kcube edges are included, to avoid
%   double counting problems.

function [klist_s, klist_r] = kmesh_gen(Ham_obj,KCUBE_BULK,kopts)
arguments
    Ham_obj {mustBeA(Ham_obj,{'HR','HK','Htrig'})}
    KCUBE_BULK double =[]; % [original_point; vk1; vk2; vk3]
    
    kopts.nk double = [1 1 1]; % [nk1 nk2 nk3]
    kopts.vk = [1 0 0; 0 1 0; 0 0 1]; % [vk1; vk2; vk3]
    kopts.original_point = [0 0 0];
    kopts.mode {mustBeMember(kopts.mode,{'corner','center'})} = 'corner'
    kopts.edge {mustBeMember(kopts.edge,{'half','full'})} = "half";
end
if ~isempty(KCUBE_BULK)
    original_point = KCUBE_BULK(1,:);
    vk = KCUBE_BULK(2:4,:);
else
    original_point = kopts.original_point;
    vk = kopts.vk;
end
nk = kopts.nk;
%%
if kopts.edge == "half"
    klist_k1 = linspace3([0 0 0], vk(1,:), nk(1)+1);
    klist_k2 = linspace3([0 0 0], vk(2,:), nk(2)+1);
    klist_k3 = linspace3([0 0 0], vk(3,:), nk(3)+1);
elseif kopts.edge == "full"
    klist_k1 = linspace3([0 0 0], vk(1,:), nk(1));
    klist_k2 = linspace3([0 0 0], vk(2,:), nk(2));
    klist_k3 = linspace3([0 0 0], vk(3,:), nk(3));    
end
%%
klist_s = zeros(nk(1) * nk(2) * nk(3), 3);
for c = 1:nk(3)
    start1 = nk(1)*nk(2) * (c-1);
    klist_s( start1+1: nk(1)*nk(2)+start1, 1:3)=...
        klist_s( start1+1: nk(1)*nk(2)+start1, 1:3) + klist_k3(c,1:3);
    for b = 1:nk(2)
        start2 = start1 + (b-1)*nk(1);
        klist_s( start2+1: nk(1)+start2, 1:3)=...
            klist_s( start2+1: nk(1)+start2, 1:3) +...
            klist_k1(1:nk(1), 1:3) + klist_k2(b, 1:3);
    end
end
klist_s = klist_s + original_point;
% dk_s = abs(cross(vk(1,:)/nk(1), vk(2,:)/nk(2)) * (vk(3,:)/nk(3))');
%%
if kopts.mode == "corner"
elseif kopts.mode == "center"
    klist_s = klist_s - (vk(1,:) + vk(2,:) + vk(3,:))./2;
end
%%
klist_r = klist_s * Ham_obj.Gk;
% vk_r = vk * Ham_obj.Gk;
% for i = 1:3
%     if nk(i) == 1
%         Id = [0 0 0];
%         Id(i) = 1;
%         vk_r(i,:) = Id;
%     end
% end
% dk_r = abs(cross(vk_r(1,:)/nk(1), vk_r(2,:)/nk(2)) * (vk_r(3,:)/nk(3))');
%%
% varargout{1} = klist_s;
% varargout{2} = klist_r;
% varargout{3} = dk_s;
% varargout{4} = dk_r;
%% 
% kstruct.klist_s = klist_s;
% kstruct.klist_r = klist_r;
% kstruct.dk_s = dk_s;
% kstruct.dk_r = dk_r;
% kstruct.nk = nk;
% varargout{5} = kstruct;
end

function vlist = linspace3(v1, v2, n)
dim = length(v1);
vlist = zeros(n, dim);
for i = 1:dim
    vlist(1:n,i) = linspace(v1(i),v2(i),n);
end
end