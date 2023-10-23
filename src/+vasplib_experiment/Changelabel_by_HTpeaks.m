function [Block_struct,fig,ax,f_handleT,f_handleH] = Changelabel_by_HTpeaks(time,T,H,Block_list,zeroH_base,refineforTandtime,headstep)
import vasplib_tool.*;
import vasplib_experiment.*
% H = PpmsVsmData.Magnetic_Field_Oe

if nargin < 5
    zeroH_base  =99;
end
if nargin < 6
    refineforTandtime = 10;
end
if nargin < 7
    headstep = 5;
end
    % save first
    time_save =  time;
    H_save = H;
    T_save = T;
    time_save = time_save - time_save(1);
    % shift
    %time(1:headstep) = [];
    H(1:headstep) = [];
    T(1:headstep) = [];
    H_refined = round(H/zeroH_base);
    %time_refined = round(time*refineforTandtime)/refineforTandtime;
    T_refined =  round(T*refineforTandtime)/refineforTandtime;   
    % ------- init
    MaxT= max(T);
    MaxH= max(H);
    nT = length(T_refined);
    
    % plot first
    [fig,ax] = subcreat_figure(2,1);
    plot(ax(1),time_save,T_save,'LineWidth',1.5,'DisplayName','T-t');
    plot(ax(2),time_save,H_save,'LineWidth',1.5,'DisplayName','H-t');
    % give xylabel
    ylabel(ax(1),'Temperature (K)');
    ylabel(ax(2),'MagField (Oe)');
    xlabel(ax(2),'Time-line (s)');
    % give title
    title(ax(1),'T-t');title(ax(2),'H-t');
    
    if nargin <4
        % ------- peaks
        % first peak
        reference_Label_list1 = reference_Label(-T_refined,'peaks');
        % second peak
        reference_Label_list2 = reference_Label(T_refined,'peaks');
        % third peak
        reference_Label_list3 = reference_Label(-H_refined,'peaks');
        % fourth peak
        reference_Label_list4 = reference_Label(H_refined,'peaks');
        % ------- gave label
        Label_list = [       1 ;...
            reference_Label_list1;...
            reference_Label_list2;...
            reference_Label_list3;...
            reference_Label_list4;
            nT-headstep];
        % Give final Data;
        Label_list = Label_list + headstep;
        Label_list = sort(Label_list);
        %disp(Label_list);
        % Give_referecnce
        stem(ax(1),time_save(Label_list),repmat(MaxT,[length(Label_list),1]));  
        stem(ax(2),time_save(Label_list),repmat(MaxH,[length(Label_list),1]));
        % make slides
        Label_block(:,1) = Label_list(1:end-1);
        Label_block(:,2) = Label_list(2:end);
        % New block 
        Block_list= Edit_Label_block(time_save,Label_block,MaxT,MaxH,ax);
    end
    % Give new block name
    Block_struct = Give_block_name(time_save,Block_list,MaxT,MaxH,ax);
    % Give final plot 
    [f_handleT,f_handleH] = New_Label_block_struct_plot(time_save,Block_struct,MaxT,MaxH,ax);
    
end






function [f_handleT,f_handleH] = New_Label_block_struct_plot(time,New_Label_block_struct,MaxT,MaxH,ax)
import vasplib_experiment.*
    axT = ax(1);axH = ax(2);
    [nLabel_block,~] = size(New_Label_block_struct);
    for i = 1:nLabel_block
        ablock = time(New_Label_block_struct(i,:).block);
        posT = [ablock(1),0,ablock(2)-ablock(1),MaxT] ;
        posH = [ablock(1),0,ablock(2)-ablock(1),MaxH] ;
        name = New_Label_block_struct(i,:).name;
        if New_Label_block_struct(i,:).type == 1
            color = 'y';
        elseif New_Label_block_struct(i,:).type == 2
            color = 'g';
        elseif New_Label_block_struct(i,:).type == 3
            color = 'r';
        else
            
        end
        [r_handleT,f_handleT(i)] = plot_rectangle(posT,axT,color,name);
        delete(r_handleT);
        [r_handleH,f_handleH(i)] = plot_rectangle(posH,axH,color,name);
        delete(r_handleH);
    end

end
function New_Label_block = Edit_Label_block(time,Label_block,MaxT,MaxH,ax)
import vasplib_experiment.*
    axT = ax(1);axH = ax(2);
    COUNT = 1;
    nLabel_block = length(Label_block);
    inputchar = 'y';
    New_Label_block=[];
    while inputchar ~= 'q' 
        if COUNT > nLabel_block 
            break;
        end
        nLabel_block = length(Label_block);
        ablock = time(Label_block(COUNT,:));
        posT=[ablock(1),0,ablock(2)-ablock(1),MaxT] ;
        posH=[ablock(1),0,ablock(2)-ablock(1),MaxH] ;
        print_information_in_Edit_Label_block(Label_block(COUNT,:),length(New_Label_block)+1);
        [r_handleT,f_handleT] = plot_rectangle(posT,axT);
        [r_handleH,f_handleH] = plot_rectangle(posH,axH);
        inputchar = input('please input your command:\n','s');
        %reminber delete handle
        if strcmp(inputchar,'y')
            New_Label_block = [New_Label_block;Label_block(COUNT,:)];
            COUNT = COUNT+1;
            delete(r_handleT);delete(r_handleH);delete(f_handleT);delete(f_handleH);
        elseif strcmp(inputchar,'n')
            COUNT = COUNT+1;
            delete(r_handleT);delete(r_handleH);delete(f_handleT);delete(f_handleH);
        elseif strcmp(inputchar,'q')
            COUNT = COUNT+1;
            delete(r_handleT);delete(r_handleH);delete(f_handleT);delete(f_handleH);
            %do nothing
            continue;
        elseif contains(inputchar,',')
            inputstring_cell = strsplit(inputchar,',');          
            inputblock(1) = double(string(inputstring_cell{1}));
            inputblock(2) = double(string(inputstring_cell{2}));
            New_Label_block =[New_Label_block;[inputblock(1) ,inputblock(2) ]];
            COUNT = COUNT+1;
            delete(r_handleT);delete(r_handleH);delete(f_handleT);delete(f_handleH);
        else
            if isempty(str2num(input_number))
                fprint('wrong input');
            else
                input_number = str2num(inputchar);
                ablock = [ablock(1),inputnumber];
                bblock = [inputnumber,ablock(2)];
                Label_block = [Label_block(1:COUNT-1,:);ablock ;bblock;Label_block(COUNT+1:end,:)];
                delete(r_handleT);delete(r_handleH);delete(f_handleT);delete(f_handleH);              
            end

            % insert
            %COUNT = COUNT;
        end
    end
end
function New_Label_block_struct = Give_block_name(time,Label_block,MaxT,MaxH,ax)
import vasplib_experiment.*
    axT = ax(1);axH = ax(2);
    COUNT = 1;
    [nLabel_block,~] = size(Label_block);
    inputchar = 'y';
     while inputchar ~= 'q' 
        if COUNT > nLabel_block 
            break;
        end
        New_Label_block_struct(COUNT,:).block = Label_block(COUNT,:);
        ablock = time(Label_block(COUNT,:));
        posT=[ablock(1),0,ablock(2)-ablock(1),MaxT] ;
        posH=[ablock(1),0,ablock(2)-ablock(1),MaxH] ;
        print_information_in_Givename_Label_block(Label_block(COUNT,:),COUNT);
        [r_handleT,f_handleT] = plot_rectangle(posT,axT);
        [r_handleH,f_handleH] = plot_rectangle(posH,axH);
        inputchar = input('','s');
        %disp(inputchar);
        if strcmp(inputchar,'q')
            delete(r_handleT);delete(r_handleH);delete(f_handleT);delete(f_handleH);
            continue;
        elseif isempty(str2num(inputchar))
            fprintf('wrong input, again');       
            continue
        elseif str2num(inputchar) == 1
%             disp('why');
            New_Label_block_struct(COUNT,:).type = 1;
        elseif str2num(inputchar) == 2
            New_Label_block_struct(COUNT,:).type = 2;           
        elseif str2num(inputchar) == 3
            New_Label_block_struct(COUNT,:).type = 3;           
        else            
        end
        inputchar2 = input('Please give a handy name for this slice:\n','s');
        New_Label_block_struct(COUNT,:).name = inputchar2;
        delete(r_handleT);delete(r_handleH);delete(f_handleT);delete(f_handleH);
        COUNT = COUNT+1;
     end  
end
function print_information_in_Edit_Label_block(a_block,iseq)
    fprintf('Here, we should check the data slice should be used or not?\n');
    fprintf('the  %d th data slice is [%d , %d]\n',iseq,a_block(1),a_block(2));
    fprintf('Use it with nochange, type: y \n');
    fprintf('Never use it        , type: n \n');
    fprintf('Quit this section   , type: q \n');
    fprintf('Modify slice        , type two new number with , (example: %d,%d) \n',...
        a_block(1)+10,a_block(2)-10);
    fprintf('Spilit slice        , type one number in the range [%d , %d], (example: %d) \n',...
        a_block(1),a_block(2),mean(a_block));
end
function print_information_in_Givename_Label_block(a_block,iseq)
    fprintf('Here, we should give the data slice property and handy name \n');
    fprintf('the  %d th data slice is [%d , %d]\n',iseq,a_block(1),a_block(2));
    fprintf('Please choose the type of experiment in this slice :\n');
    fprintf('FC, type: 1 \n');
    fprintf('ZFC, type: 2 \n');
    fprintf('MH, type: 3 \n');
    fprintf('Quit this section   , type: q \n');
end
function reference_Label_list = reference_Label(data,mode)
    if strcmp(mode,'peaks')
        [~,reference_Label_list] = findpeaks(smooth( data));
    end
end


%     %%
% 
%     %
%     for i = 1:nT_dif_t
%         temp_T_diff = T_dif_t_floor(i);
%         if abs(temp_T_diff) <tinysmall
%             if seq_label >0
%                 tinysmall_COUNTCOUNT = tinysmall_COUNTCOUNT+1;
%             end
%             seq_label = seq_label+1;
%             change_label = 0;
%             temp_T_diff = 0;
%         else
%             if change_label >0
%                 Change_COUNTCOUNT = Change_COUNTCOUNT+1;
%             end
%             seq_label = 0;
%             change_label = change_label + 1;
%         end
%         
%         if temp_T_diff*signlabel <0
%             % change sign
%             signlabel = signlabel*-1;            
%             COUNTCOUNT = COUNTCOUNT + 1 ;
%             Label_list(COUNTCOUNT) = i;
%         elseif tinysmall_COUNTCOUNT > NumZeroT & Allowconutlabel == 1
%             COUNTCOUNT = COUNTCOUNT + 1 ;
%             Label_list(COUNTCOUNT) = i-tinysmall_COUNTCOUNT;
%             tinysmall_COUNTCOUNT = 0;
%             Allowconutlabel = 0;
%             % no signchange, do nothing     
%         end
%         
%         if Change_COUNTCOUNT > NumZeroT 
%             Allowconutlabel = 1;
%             Change_COUNTCOUNT = 0;
%         end
%     end

    %%  second test H
%     MaxH= max(H);
%     tinysmall_COUNTCOUNT = 0;
%     Change_COUNTCOUNT = 0;
%     signlabel = 1 ;
%     seq_label = 1;
%     change_label = 1;
%     nH_dif_t = length(H_dif_t);
%     Allowconutlabel = 1;
%     % nH = nH_dif_t +1
%     nH = length(H_refined);
%     %%
% 
%     %
%     for i = 1:nH_dif_t
%         temp_H_diff = H_dif_t_floor(i);
%         temp_H = H_refined(i);
%         if abs(temp_H_diff) <tinysmall
%             if seq_label >0
%                 tinysmall_COUNTCOUNT = tinysmall_COUNTCOUNT+1;
%             end
%             seq_label = seq_label+1;
%             change_label = 0;
%             temp_H_diff = 0;
%         else
%             if change_label >0
%                 Change_COUNTCOUNT = Change_COUNTCOUNT+1;
%             end
%             seq_label = 0;
%             change_label = change_label + 1;
%         end
%         
%         if temp_H_diff*signlabel <0
%             % change sign
%             signlabel = signlabel*-1;            
%             COUNTCOUNT = COUNTCOUNT + 1 ;
%             Label_list(COUNTCOUNT) = i;
%         elseif tinysmall_COUNTCOUNT > NumZeroT & Allowconutlabel == 1
%             COUNTCOUNT = COUNTCOUNT + 1 ;
%             Label_list(COUNTCOUNT) = i-tinysmall_COUNTCOUNT;
%             tinysmall_COUNTCOUNT = 0;
%             Allowconutlabel = 0;
%         else
%             % no signchange, do nothing     
%         end
%         if temp_H == 0
%             COUNTCOUNT = COUNTCOUNT + 1 ;
%             Label_list(COUNTCOUNT) = i;
%         end
%         if Change_COUNTCOUNT > NumZeroT 
%             Allowconutlabel = 1;
%             Change_COUNTCOUNT = 0;
%         end
%     end