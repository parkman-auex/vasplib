

load bulkek_plane-matlab.dat

style = 2;
alpha = 0.7;

Kx=bulkek_plane_matlab(:,1);
Ky=bulkek_plane_matlab(:,2);
Kz=bulkek_plane_matlab(:,3);
K1=bulkek_plane_matlab(:,4);
K2=bulkek_plane_matlab(:,5);
K3=bulkek_plane_matlab(:,6);

size =sqrt(length(Kx));

Eoccu_m1 = reshape(bulkek_plane_matlab(:,7),size,size);
Eoccu = reshape(bulkek_plane_matlab(:,8),size,size);
Eoccu_p1 = reshape(bulkek_plane_matlab(:,9),size,size);
Eoccu_p2 = reshape(bulkek_plane_matlab(:,10),size,size);

KXX = reshape(Kx,size,size);
KYY = reshape(Ky,size,size);
KZZ = reshape(Kz,size,size);
K11 = reshape(K1,size,size);
K22 = reshape(K2,size,size);
K33 = reshape(K3,size,size);
if style ==1
surf(KXX,KYY,Eoccu_m1,'linestyle','none','FaceAlpha',alpha);
hold on;
surf(KXX,KYY,Eoccu,'linestyle','none','FaceAlpha',alpha);
hold on;
surf(KXX,KYY,Eoccu_p1,'linestyle','none','FaceAlpha',alpha);
hold on;
surf(KXX,KYY,Eoccu_p2,'linestyle','none','FaceAlpha',alpha);
hold on;
elseif style ==2
surf(K11,K22,Eoccu_m1,'linestyle','none','FaceAlpha',alpha);
hold on;
surf(K11,K22,Eoccu,'linestyle','none','FaceAlpha',alpha);
hold on;
surf(K11,K22,Eoccu_p1,'linestyle','none','FaceAlpha',alpha);
hold on;
surf(K11,K22,Eoccu_p2,'linestyle','none','FaceAlpha',alpha);
hold on;      
end
