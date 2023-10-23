function ax = ShowSurfIntegral(Gk,klist_r,dk1L,dk2L,dSumL,options)
arguments
    Gk =[];
    klist_r =[];
    dk1L =[];
    dk2L =[];
    dSumL =[];
    options.dName = 'dBerryphase';
    options.SumName = 'Berryphase';
    options.FinalName = 'ChernNumber = ';
    options.oneshot = true;
    options.view = [0,90];
end
fig = figure('PaperType','a4letter','PaperSize',[16 8],'Color','white','Units','normalized','Position',[0.1,0.1,0.8,0.6]);
ax0 = subplot(1,3,1,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
hold(ax0,'all');
[fig,ax0] = vasplib.BZplot(Gk,'Gk',true,'ax',ax0,'fig',fig,'color','r','alpha',0.3);
title(ax0,'integral patch');
axis(ax0,'equal');
view(ax0,options.view(1),options.view(2));
ax1 = subplot(1,3,2,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
hold(ax1,'all');
xlabel(ax1,'kpath num');
ylabel(ax1,options.dName);
ax2 = subplot(1,3,3,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
hold(ax2,'all');
xlabel(ax2,'kpath num');
ylabel(ax2,options.SumName);
kn = size(klist_r,1);
%             if length(dk1L) == 1
%                 dk1L = repmat(dk1L,[kn,1]);
%             end
%             if length(dk2L) == 1
%                 dk2L = repmat(dk2L,[kn,1]);
%             end
dk1dk2L = dk1L+dk2L;
KsmallL = zeros(kn,3,4);
KsmallL(:,:,1) = klist_r;
KsmallL(:,:,2) = klist_r+dk1L;
KsmallL(:,:,3) = klist_r+dk2L;
KsmallL(:,:,4) = klist_r+dk1dk2L;
SUM = 0;
dSUM = 0;
dSUM_old = dSUM;
count = 0;
dSumL = real(dSumL);
for ki = 1:kn
    K = klist_r(ki,:);
    dSUM =dSumL(ki);
    count = count +1;
    Ksmall = [KsmallL(ki,:,1);KsmallL(ki,:,2);KsmallL(ki,:,4);KsmallL(ki,:,3)];
    patch(ax0,Ksmall(:,1),Ksmall(:,2),Ksmall(:,3),'g','EdgeColor','none');
    %
    xlabel(ax1,string(K(1))+","+string(K(2))+","+string(K(3)));
    area(ax1,[count-1,count],[dSUM_old,dSUM]);
    dSUM_old= dSUM;
    %
    SUM = SUM+dSUM;
    stem(ax2,count,SUM);
    if options.oneshot
    else
        drawnow;
    end
end
title(ax2,options.FinalName+string(SUM));
axis(ax0,'equal');
ax = [ax0,ax1,ax2];
end
