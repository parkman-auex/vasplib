
function [num_bands,Nmin,Nmax] = selectbands(EIGENCAR,Rm,mode)
    if nargin <3
        mode = 'win_gen';
    end
    disp("please input Nmin and Nmax :");
    Nmin=input('please input the Nmin\n');
    Nmax=input('please input the Nmax\n');
    if strcmp(mode,'win_gen')
        num_bands=Nmax-Nmin+1;
        [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm);
        size_=[0.2 0.2 0.6 0.6];
        fontname="Helvetica";
        %%
        kpathname=kpoints_name;
        %% listK
        klist=klist_l;
        %% paper set
        figure_k = figure('PaperType','a4letter','PaperSize',[8 6],'Color',[1 1 1]);
        axes1 = axes('Parent',figure_k ,'Position',size_,'LineWidth',1.5,'FontSize',30.0,'FontName',fontname);         %Position [left bottom width height]
        box(axes1,'on');
        hold(axes1,'all');
        %% plot whole
        Nbands=size(EIGENCAR,1);
        for Ei=1:Nbands
            plot(klist,EIGENCAR(Ei,:),'LineWidth',2.0,'Color',[0 0 1],'DisplayName',num2str(Ei));
            hold on
        end
        %% plot special
        plot(klist,EIGENCAR(Nmin,:),'LineWidth',3.0,'Color',[1 0 0],'DisplayName',"Nmin");
        hold on
        plot(klist,EIGENCAR(Nmax,:),'LineWidth',3.0,'Color',[1 0 0],'DisplayName',"Nmax");
        hold on
        plot(klist,EIGENCAR(Nmin-1,:),'LineWidth',3.0,'Color',[1 1 0],'DisplayName',"Nmin",'Linestyle',"--");
        hold on
        plot(klist,EIGENCAR(Nmax+1,:),'LineWidth',3.0,'Color',[1 1 0],'DisplayName',"Nmax",'Linestyle',"--");
        hold on
        %% plot set
        xmin=kpoints_l(1);                         %\
        xmax=kpoints_l(length(kpoints_l));
        ymin=min(min(EIGENCAR(:,:)));
        ymax=max(max(EIGENCAR(:,:)));
        axis([xmin xmax ymin ymax]);
        % pre set pathname
        for i=1:length(kpathname)
            kpathname(i)=strrep(kpathname(i),'GAMMA','G');
            kpathname(i)=strrep(kpathname(i),'G','\Gamma');
        end
        %% set label
        set(gca,'XTickLabel','')  ;
        set(axes1,'FontName','Helvetica','FontSize',30,'LineWidth',1.5,...
            'XTick',kpoints_l,'XTickLabel',...
            kpathname);
        % Create ylabel
        ylabel('E-E_F (eV)','FontWeight','bold','FontName',fontname);
        % regrence line
        X=[xmin xmax];
        Y=[ymin ymax];
        line(X,[0 0],'LineWidth',0.5,'Color',[0 0 0],'DisplayName','fermi-level')
        for i=1:1: (length(kpathname)-2)
            X=[kpoints_l(i+1) kpoints_l(i+1)];
            Y=[ymin ymax];
            line(X,Y,'LineWidth',0.2,'Color',[0 0 0],'DisplayName','K-path');
        end
        title("check Erange");
    elseif strcmp(mode,'dis_check')
        infinite_small = 1e-2;
        dis_min = min(EIGENCAR(Nmin,:))-infinite_small;
        dis_max = max(EIGENCAR(Nmax,:))+infinite_small;
        [allk_min,anyk_max]=findbands(EIGENCAR,[dis_min,dis_max]);
        foz_min = max(EIGENCAR(Nmin,:))+infinite_small;
        foz_max = min(EIGENCAR(Nmax,:))-infinite_small;
        disp('recomand:')
        fprintf('Nbands = %d !',anyk_max.nbands);
        fprintf('choose (%d - %d)\n',anyk_max.n1,anyk_max.n2);
        fprintf('dis_win_min = %f\n',dis_min);
        fprintf('dis_win_max = %f\n',dis_max);
        fprintf('dis_froz_min = %f\n',foz_min);
        fprintf('dis_froz_max = %f\n',foz_max);
        fprintf('another choice: \n')
        fprintf('dis_win_min < %f\n',dis_min);
        fprintf('dis_win_max > %f\n',dis_max);
        fprintf('dis_froz_min > %f\n',foz_min);
        fprintf('dis_froz_max < %f\n',foz_max);        
    end
end
