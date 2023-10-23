function kpath_card_gen()

%% KPOINT_path card 


for flag=1:1
%note : ['Line-Mode' ],['Reciprocal' ]}
 % get whole information
    KPOINTS=fopen('KPOINTS','r');                 %open KPOINTS
    temp_i=1;
    while ~feof(KPOINTS)
        KPOINTS_information{temp_i}=fgets(KPOINTS);
        temp_i=temp_i+1;
    end
    
    fclose(KPOINTS);                              %close KPOINTS
    % nodes
    nodes = str2num(KPOINTS_information{2});
    % kpoints list
    headline=5;
    kpoints=[];
    kpoints_name_pre=[];
    for i=headline:temp_i-1
        kpoints_line=strtrim(KPOINTS_information{i});
        kpoints_line=regexp(kpoints_line, '\s+', 'split');
         if length(kpoints_line)>1      
             kpoints=[kpoints;[str2num(kpoints_line{1}),str2num(kpoints_line{2}),str2num(kpoints_line{3})]];
             kpoints_name_pre=[kpoints_name_pre;string(kpoints_line{4})];% 
         end
    end
% wannier90 k=path card ____________________________________________________________________________________________________________
    fid = fopen('wannier90kpath_card','w');
    % function string
    head_string="begin kpoint_path ";
    end_string ="end kpoint_path   ";
    fprintf(fid,"%s\n",head_string);
    for i=1:length(kpoints)/2
        fprintf(fid,"%1s ",strrep(kpoints_name_pre(2*i-1),'GAMMA','G'));
        fprintf(fid,"%9s ",num2str(kpoints(2*i-1,1)));
        fprintf(fid,"%9s ",num2str(kpoints(2*i-1,2)));
        fprintf(fid,"%9s ",num2str(kpoints(2*i-1,3)));
        fprintf(fid,"%1s ",strrep(kpoints_name_pre(2*i),'GAMMA','G'));
        fprintf(fid,"%9s ",num2str(kpoints(2*i,1)));
        fprintf(fid,"%9s ",num2str(kpoints(2*i,2)));
        fprintf(fid,"%9s \n",num2str(kpoints(2*i,3)));
    end
    fprintf(fid,"%s\n",end_string);
    fclose(fid);
end
% wt.in k=path card ____________________________________________________________________________________________________________
    fid = fopen('wt_in_kpath_card','w');
    % function string
    head_string="KPATH_BULK            ! k point path ";
    fprintf(fid,"%s\n",head_string);
    fprintf(fid,"%d\n",length(kpoints)/2); 
    for i=1:length(kpoints)/2
        fprintf(fid,"%1s ",strrep(kpoints_name_pre(2*i-1),'GAMMA','G'));
        fprintf(fid,"%9s ",num2str(kpoints(2*i-1,1)));
        fprintf(fid,"%9s ",num2str(kpoints(2*i-1,2)));
        fprintf(fid,"%9s ",num2str(kpoints(2*i-1,3)));
        fprintf(fid,"%1s ",strrep(kpoints_name_pre(2*i),'GAMMA','G'));
        fprintf(fid,"%9s ",num2str(kpoints(2*i,1)));
        fprintf(fid,"%9s ",num2str(kpoints(2*i,2)));
        fprintf(fid,"%9s \n",num2str(kpoints(2*i,3)));
    end

    fclose(fid);
end
