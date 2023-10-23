%%  DOS_plot
%
%  creat a list of DOS
%
%
% * Label:
%
%%
%% Usage:
%
% * [figure_dos,Pdensity,Totdensity] = dosplot(Ecut,choose_list,mode)
% * [figure_dos,Pdensity,Totdensity] = dosplot(Ecut,choose_list)
% * [figure_dos,Pdensity,Totdensity] = dosplot(Ecut)
%
%% Input:
%
% # Ecut:
% # choose_list:
% # mode:
%
%% Output:
%
% # figure_dos:
% # Pdensity:
% # Totdensity:
%
%% example:
%   commmad
%
%
%   result
%
%% Note:
%
%  This function is not effienct , need to be improved!!!.
%
%% Change log
%
% * Document Date: 2020/12/03
% * Creation Date: 2020/12/03
% * Last updated : 2022/02/23
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code :
%
function [figure_dos,Pdensity,Totdensity] = dosplot(Ecut,choose_list,mode,options)
    arguments
        Ecut = [-2,2];
        choose_list = -1;
        mode = 'vaspkit-dir';
        options.d = true;
    end
    %--------  init  --------
    
    import park.*
    import vasplib_tool.*
    %--------  narg  --------
    %--------  chek  --------
    %--------  juge  --------
    if strcmp(mode,'vaspkit-dir')
        % dat gen TDOS
        filename = 'tdos.dat';
        [Totdensity,EIGENCAR,~] = WEIGHTCAR_gen(filename,-1,'vaspkit-DOS-silence');
        % dat gen PDOS
        [Pdensity,Pdos_namelist,width] = WEIGHTCAR_gen_dir('DOS');
        %[~,width] = size(Pdensity(:,:,1));
        % f or not
        if width >10
            f_mode = true;
            num2orbital_name =  orbital_maprule_gen(1);
        else
            f_mode = false;
            num2orbital_name = orbital_maprule_gen(0);
        end
        if choose_list <0
            %
            fprintf('vaspkit DOS dir \n');
            % Tot first
            ax(1) = dos_plot_set(EIGENCAR,Totdensity,'DOS-TOT',Ecut,'DOS-TOT');
            figure_dos(1) = ax(1).Parent;
            % each element second
            EIGENCAR_element(:,1)= EIGENCAR;
            WEIGHTCAR_element(:,1)= Totdensity;
            Name_list(1,:)="TOT";
            for i =1:length(Pdos_namelist)
                EIGENCAR_element(:,i+1) = EIGENCAR;
                WEIGHTCAR_element(:,i+1) = Pdensity(:,end,i);
                Name_list(i+1,:) = Pdos_namelist(i,:);
            end
            [ax(2)] = dos_plot_set( EIGENCAR_element,WEIGHTCAR_element,Name_list,Ecut,'DOS-element');figure_dos(2)= ax(2).Parent;
            % each orbital
            Name_list = [""];
            EIGENCAR_orbital(:,1)= EIGENCAR;
            WEIGHTCAR_orbital(:,1)= Totdensity;
            Name_list(1,:)='TOT';
            % s
            EIGENCAR_orbital(:,2)= EIGENCAR;
            temp_WEIGHTCAR = [];
            for i =1:length(Pdos_namelist)
                temp_WEIGHTCAR = [ temp_WEIGHTCAR,Pdensity(:,1,i)];
            end
            WEIGHTCAR_orbital(:,2)= sum(temp_WEIGHTCAR,2);
            Name_list(2,:)='s';
            % p
            EIGENCAR_orbital(:,3)= EIGENCAR;
            temp_WEIGHTCAR = [];
            for i =1:length(Pdos_namelist)
                temp_WEIGHTCAR = [ temp_WEIGHTCAR,Pdensity(:,2:4,i)];
            end
            WEIGHTCAR_orbital(:,3)= sum(temp_WEIGHTCAR,2);
            Name_list(3,:)='p';
            % d
            if options.d
                EIGENCAR_orbital(:,4)= EIGENCAR;
                temp_WEIGHTCAR = [];
                for i =1:length(Pdos_namelist)
                    temp_WEIGHTCAR = [ temp_WEIGHTCAR,Pdensity(:,5:9,i)];
                end
                WEIGHTCAR_orbital(:,4)= sum(temp_WEIGHTCAR,2);
                Name_list(4,:)='d';
            end
            % f
            if width >10
                EIGENCAR_orbital(:,5)= EIGENCAR;
                temp_WEIGHTCAR = [];
                for i =1:length(Pdos_namelist)
                    temp_WEIGHTCAR = [ temp_WEIGHTCAR,Pdensity(:,10:16,i)];
                end
                WEIGHTCAR_orbital(:,5)= sum(temp_WEIGHTCAR,2);
                Name_list(5,:)='f';
            end
            [ax(3)] = dos_plot_set( EIGENCAR_orbital,WEIGHTCAR_orbital,Name_list,Ecut,'DOS-orbital');figure_dos(3)= ax(3).Parent;
            % each element
            count = 4;
            if options.d
            else
                width = 4;
            end
            for i = 1:length(Pdos_namelist)
                Name_list = [""];
                EIGENCAR_atom(:,1)= EIGENCAR;
                WEIGHTCAR_atom(:,1)= Totdensity;
                Name_list(1,:) = "TOT"';
                WEIGHTCAR_atom(:,2:width+1) = Pdensity(:,1:width,i);
                for j =1:width
                    EIGENCAR_atom(:,j+1) = EIGENCAR;

                    tempstring = strcat(strcat(Pdos_namelist(i,:),'-'),num2orbital_name(j));
                    %                     disp(num2orbital_name(j));
                    %                     disp(tempstring);
                    Name_list(j+1,:) = tempstring;
                end
                figure_name = strcat('DOS-orbital-',Pdos_namelist(i,:));
                [ax(count)] = dos_plot_set(EIGENCAR_atom,...
                    WEIGHTCAR_atom,Name_list,Ecut, figure_name);figure_dos(count) = ax(count).Parent;
                count = count+1;
            end
            % Your choose
        else
            % single mode
            print_PDOS_information_in_dir(Pdos_namelist,width);
            if ISA(choose_list,'double')
                fprintf(" Let's do it maully ");
            else
                EIGENCAR_atom(:,1)= EIGENCAR;
                WEIGHTCAR_atom(:,1)= Totdensity;
                Name_list(1,:) = "TOT"';
                Nchoose = length(choose_list);
                for i = 1:Nchoose
                    EIGENCAR_atom(:,i+1) = EIGENCAR;
                    Name_list(i+1,:) = choose_list(i).displayname;
                    WEIGHTCAR_atom(:,i+1) = sum(sum(Pdensity(:,...
                        choose_list(i).obital_seq_list,...
                        choose_list(i).atom_seq_list),3),2);
                end
                figure_name = 'PDOS';
                [ax] = dos_plot_set(EIGENCAR_atom,...
                    WEIGHTCAR_atom,Name_list,Ecut, figure_name);
                figure_dos = ax.Parent;
            end
        end
    end


end




% function [figure_k,axis] = dos_plot_one(EIGENCAR,WEIGHTCAR,Name,axis,figure_k)
%     Fontname="Helvetica";
%     if nargin < 4
%         figure_k= figure('PaperType','a4letter','PaperSize',[8 8],'Color','white');
%         axis = axes('Parent',figure_k,'LineWidth',1.5,'FontSize',24.0,'FontName',Fontname);         %Position [left bottom width height]
%     end
%     plot(axis,EIGENCAR,WEIGHTCAR,'LineWidth',2.0,'color',[rand,rand,rand],'DisplayName',Name);
% end