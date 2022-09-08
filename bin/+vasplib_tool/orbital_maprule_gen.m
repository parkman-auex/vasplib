function maprule= orbital_maprule_gen(f_or_not)
if nargin <1
    f_or_not = 0;
end
if f_or_not == 1
    maprule= containers.Map({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17},...
        {'s','py','pz','px','dxy','dyz','dz2','dxz','dx2-y2','fy3x2','fxyz','fyz2','fz3','fxz2','fzx2','fx3','tot'});
elseif f_or_not  == 0
    maprule= containers.Map({1,2,3,4,5,6,7,8,9,10},...
        {'s','py','pz','px','dxy','dyz','dz2','dxz','dx2-y2','tot'});
end
end