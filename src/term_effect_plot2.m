function [klist]=term_effect_plot2(EIGENCAR,termlist,termname,Ecut,size_,fontname,titlestring)
%% set default para
if nargin ==1
    termlist=linspace(-1.5,1.5,100);
    Ecut=[-10,10];
    size_=[0.2 0.2 0.6 0.6];
    fontname="Times New Roman";
    termname="tb/ta";
    
elseif nargin==2
    termname="ta/tb";
    Ecut=[-10,10];
    size_=[0.2 0.2 0.6 0.6];
    fontname="Times New Roman";

elseif nargin==3        
    Ecut=[-2,2];
    size_=[0.2 0.2 0.6 0.6];
    fontname="Times New Roman";
elseif nargin==4        
    size_=[0.2 0.2 0.6 0.6];
    fontname="Times New Roman";
elseif nargin==5        
    fontname="Times New Roman";
end
if nargin<6
    % Create title
    dirname=pwd;
    dirname=strsplit(dirname,'/');
    titlestring=dirname(length(dirname));
end

%% listK

format long;


nterm=length(termlist);
klist=[];

%% paper set 
A4 = figure('PaperType','a4letter','PaperSize',[8 6],'Color',[1 1 1]);
axes1 = axes('Parent',A4,'Position',size_,'LineWidth',1.5,'FontSize',30.0,'FontName',fontname);         %Position [left bottom width height]
box(axes1,'on');
hold(axes1,'all');
%% plot whole

    Nbands=size((EIGENCAR),1);
    for En=1:Nbands
        plot(termlist,EIGENCAR(En,:),'LineWidth',1.0,'Color',[0 0 0],'DisplayName',num2str(En)+"bands",'Color',[0 0 1]);
    hold on
    end

%% plot set 
    xmin=termlist(1);                         %\
    xmax=termlist(nterm);
    ymin=Ecut(1);
    ymax=Ecut(2);
    axis([xmin xmax ymin ymax]); 
%% set label
% Create xlabel
    xlabel(termname,'FontWeight','bold','FontName',fontname);
    h=termlist(1);
% 
%     set(gca,'XTickLabel','')  ;                                                      
    xlabel(termname,'FontWeight','bold','FontName',fontname);
%     set(axes1,'FontName','Kozuka Gothic Pr6N','FontSize',30,'LineWidth',1.5,...
%     'XTick',[h h/2 0 -h/2 -h],'XTickLabel',...
%     ["-1" "-1/2" "0" "1/2" "1"]);
% Create ylabel
    ylabel('E(eV)','FontWeight','bold','FontName',fontname);

% regrence line    
    X=[xmin xmax];                                                
%     Y=[ymin ymax];
% line(X,[0 0],'LineWidth',0.5,'Color',[0 0 0],'DisplayName','fermi-level')
%     for i=1:1: (length(kpathname)-2)                                           
%         X=[KPOINTS_t(i+1) KPOINTS_t(i+1)];
%         Y=[ymin ymax];
%         line(X,Y,'LineWidth',0.2,'Color',[0 0 0],'DisplayName','k-point');
%     end

    title(titlestring);

end
