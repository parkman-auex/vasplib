%% main-fourband model




function [EIGENCAR_temp,WAVECAR]=ta_tb_effect(H_xyz,internum,mode,t_c)
%%init


    
if nargin<3   
    mode = 'TB' ;
end    
if nargin<4    
    t_c =0 ;
end  
if nargin<2    
    internum  = 100;
end   
    t_d =0 ;
    t_e =0 ;
%% term 
    t_b =1;
    % gen ta
    t_a_list =linspace(-2*t_b,2*t_b,internum);
    ta_n=length(t_a_list);
if strcmp(mode,'TB')
    POSCAR_read;
    [~,klist_l,klist_s,kpoints_l]=kpathgen3D(Rm);
    k_n=length(klist_l);
for i=1:ta_n
    t_a = t_a_list(i);
    %EIGENCAR_temp=[];  
    save('parm','t_a','t_b');
    H_temp = Subsall(H_xyz,'file');
    [EIGENCAR_temp{i},~] = EIGENCAR_gen(H_temp ,'t',klist_s);
    WAVECAR = [];
    %[figure_k,axes1]=bandplot(EIGENCAR,Ecut,titlestring,klist_l,kpoints_l,kpoints_name,size_,fontname,figure_k,axes1)
end



elseif strcmp(mode,'TB_fin')
    for i=1:ta_n
        t_a = t_a_list(i);
        %EIGENCAR_temp=[];
        save('parm','t_a','t_b');
        H_temp = Subsall(H_xyz,'file');
        [EIGENCAR_temp(:,i),WAVECAR{i}] = EIGENCAR_gen(H_temp ,'m',[0 0 0]);
    end
    [Nbands,~] = size(EIGENCAR_temp);
    fontname="Kozuka Gothic Pr6N";
    figure_k = figure('PaperType','a4letter','PaperSize',[8 6],'Color','white');
    axes1 = axes('Parent',figure_k ,'LineWidth',1.5,'FontSize',24.0,'FontName',fontname); 
        box(axes1,'on');
    hold(axes1,'all');
    for Ei=1:Nbands
        plot(t_a_list,EIGENCAR_temp(Ei,:),'LineWidth',1.0,'Color',[0 0 1],'DisplayName',num2str(Ei));
        hold on
    end
end


end
%%[klist]=bandplot(EIGENCAR,KPOINTS_s,Kpathname,Ecut,size_,nodes,fontname,titlestring);
%term_effect_plot(EIGENCAR,ta)
%[klist]=bandplot(EIGENCAR{nsite},KPOINTS_s,Kpathname,Ecut,size_,nodes,"Kozuka Gothic Pr6N","TB_{bulk band}(ta/tb=)");