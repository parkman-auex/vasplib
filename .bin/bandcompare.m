function figure_k=bandcompare(EIGENCAR,EIGENCAR2,Ecut,titlestring,klist_l,kpoints_l,kpoints_name,size_,fontname)
%% General
% usage figure_k=bandplot(EIGENCAR,Ecut,klist_l,titlestring,kpoints_l,kpoints_name,size_,fontname)

bandflag=1;
size_=[0.2 0.2 0.6 0.6];

fontname="Kozuka Gothic Pr6N";  
if nargin <1
    EIGENCAR=EIGENVAl_read();
    Ecut=[-2,2];
    POSCAR_read;
    [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm);   
    size_=[0.2 0.2 0.6 0.6];
    fontname="Kozuka Gothic Pr6N";  
end

%% set default para
if nargin ==2
    Ecut=[-3,3];
    POSCAR_read;
    [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm);   
    size_=[0.2 0.2 0.6 0.6];
    fontname="Kozuka Gothic Pr6N";    
elseif nargin==3
    POSCAR_read;
    [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm);   
    size_=[0.2 0.2 0.6 0.6];
    fontname="Kozuka Gothic Pr6N";
elseif nargin==5
    bandflag=0;
elseif nargin==6        
    size_=[0.2 0.2 0.6 0.6];
    fontname="Kozuka Gothic Pr6N";
    kpathnumber=length(kpoints_l);
    for j=1:kpathnumber
        tempstr="HighK"+num2str(j);
        kpathname=[kpathname tempstr];
    end
    kpoints_name=kpathname;
elseif nargin==7        
    size_=[0.2 0.2 0.6 0.6];
    fontname="Kozuka Gothic Pr6N";
elseif nargin==8        
    fontname="Kozuka Gothic Pr6N";
end
if nargin<4
    % Create title
    dirname=pwd;
    dirname=strsplit(dirname,'/');
    titlestring=dirname(length(dirname));
end
if nargin<5
    POSCAR_read;
    [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm);   
    size_=[0.2 0.2 0.6 0.6];
    fontname="Kozuka Gothic Pr6N"; 
end
%%
 %  kpathname=;
%% listK
   klist=klist_l;
%% paper set 
    figure_k = figure('PaperType','a4letter','PaperSize',[8 6],'Color',[1 1 1]);
    axes1 = axes('Parent',figure_k ,'Position',size_,'LineWidth',1.5,'FontSize',24.0,'FontName',fontname);         %Position [left bottom width height]
    box(axes1,'on');
    hold(axes1,'all');
%% plot whole
    Nbands=size(EIGENCAR,1);
    for Ei=1:Nbands
        plot(klist,EIGENCAR(Ei,:),'LineWidth',1.0,'Color',[0 0 1],'DisplayName',num2str(Ei));
        hold on
    end
    Nbands2=size(EIGENCAR2,1);
    for Ei=1:Nbands2
        plot(klist,EIGENCAR2(Ei,:),'.','Color',[1 0 0],'DisplayName',num2str(Ei)+"_Second");
        hold on
    end
            
    xmin=klist_l(1);                         %\
    xmax=klist_l(length(klist_l));
    ymin=Ecut(1);
    ymax=Ecut(2);
    axis([xmin xmax ymin ymax]);   
if bandflag ==1
%% plot set 

    
% pre set pathname
  for i=1:length(kpoints_name)
      kpoints_name(i)=strrep(kpoints_name(i),'GAMMA','G');
      %kpathname(i)=strrep(kpathname(i),'G','\Gamma');
  end
%% set label
    set(gca,'XTickLabel','')  ;
%     set(gca,'color','none');
%     set(gcf,'color','none');
%     set(gcf,'InvertHardCopy','off');
    set(axes1,'FontName','Kozuka Gothic Pr6N','FontSize',24,'LineWidth',1,...
    'XTick',kpoints_l,'XTickLabel',...
    kpoints_name);
% Create ylabel
    ylabel('E-E_F (eV)','FontWeight','bold','FontName',fontname);

% regrence line    
    X=[xmin xmax];                                                
    Y=[ymin ymax];
    line(X,[0 0],'LineWidth',0.5,'Color',[0 0 0],'DisplayName','fermi-level')
    for i=1:1: (length(kpoints_name)-2)                                           
        X=[kpoints_l(i+1) kpoints_l(i+1)];
        Y=[ymin ymax];
        line(X,Y,'LineWidth',0.2,'Color',[0 0 0],'DisplayName','K-path');
    end

%     title(titlestring);
%     titlestring2=titlestring+".eps";
%     saveas(figure_k,titlestring2);
%     tempstring="print -depsc2 "+titlestring2;
%     eval(tempstring); 
%     set(figure_k,'color','white')
%     set(gca,'color','white');
%     set(gcf,'color','white');
   
end


end
