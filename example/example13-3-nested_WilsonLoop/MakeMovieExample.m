v = VideoWriter('k_z_changingWilsonLoop-Hex.mp4','MPEG-4');
open(v);
for i= 1:181
% for i= 1:91
k_z = (i-1)/360;
%[BFCAR,~,klist_l] = BHZ_TB_n.WilsonLoop('kstart',[-1,0,k_z],'kevolution',[1,0,0]);
[BFCAR,~,klist_l] = BHZ_TB_n.WilsonLoop('kstart',[0,0,k_z]);
[fig,ax] = vasplib.WilsonLoopPlot(BFCAR,klist_l,'Color','r');
title (ax,['WilsonLoop k_z = (',num2str(i-1),'/180) \pi']); 
%   fig.Position = [391,432,758,701];
    writeVideo(v,getframe(fig));
    close(fig);%关闭句柄为H的figure；
end
close(v);