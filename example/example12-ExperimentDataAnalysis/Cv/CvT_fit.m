function value = CvT_fit(para,Debyecoe)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
powerlaw_coe = para.powerlaw_coe ;
powerlaw_index= para.powerlaw_index;
Nulear_coe  = para.Nulear_coe ;
if Debyecoe == 0
Debye_coe = para.Debye_coe ;
else
Debye_coe = Debyecoe ;
end
Cv = evalin('base','Cv');
T  = evalin('base','T');
% Cv = Cv*10^3;
Cv_prime = Nulear_coe*(T.^2)+powerlaw_coe*T.^(powerlaw_index)+Debye_coe *T.^3;

value = sum((Cv_prime-Cv).^2);
end

