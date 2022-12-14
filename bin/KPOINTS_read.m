%% KPOINTS_read
%
% get kpoints information from vasp and others
% * Label:
%
%% Description of the Function:
%%
%% Usage: 
%
% * [kpoints,nodes,kpoints_name] = KPOINTS_read(filename,mode)
% * [kpoints,nodes,kpoints_name] = KPOINTS_read(filename)
% * [kpoints,nodes,kpoints_name] = KPOINTS_read()
%
%% Input:
%  
% # input1:
% # input2:
% # input3:
%
%% Output:
%
% # output1:
% # output2:
% # output3:
%
%% example:
%   commmad
% 
% 
%   result
%   
%% Note: 
%
%  Take advantage of the scope of application of the function.
%
%% Change log
%
% * Document Date: 2020/12/03
% * Creation Date: 2020/12/03
% * Last updated : 2020/12/03
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code : 
%
function [kpoints,nodes,kpoints_name] = KPOINTS_read(filename,mode)
%--------  init  --------
%--------  narg  --------
    if nargin < 2
        mode = 'vasp';
    end
    if nargin < 1
        filename = 'KPOINTS';
    end
%--------  chek  --------
    if ~exist(filename,'file')
        fprintf('No such file: %s\n',filename)
        error();
    end
%--------  fbug  --------

%--------  juge  --------
    if strcmp(mode,'vasp')
        % get whole information
        KPOINTS=fopen(filename,'r');                 %open KPOINTS
        temp_i=0;
        while ~feof(KPOINTS)
            temp_i=temp_i+1;
            KPOINTS_information{temp_i}=fgets(KPOINTS);
        end
        fclose(KPOINTS);                              %close KPOINTS
        % nodes
        nodes = str2num(KPOINTS_information{2});
        % kpoints list
        headline=5;
        kpoints=[];
        kpoints_name_pre=[];
        for i=headline:temp_i
            kpoints_line=strtrim(KPOINTS_information{i});
            kpoints_line=regexp(kpoints_line, '\s+', 'split');
            if length(kpoints_line)>2
                kpoints=[kpoints;...
                    [str2double(kpoints_line{1}),str2double(kpoints_line{2}),str2double(kpoints_line{3})]];
                if length(kpoints_line)>3
                    kpoints_name_pre=[kpoints_name_pre;...
                        string(kpoints_line{4})];%
                else
                    % we consider to make a BZ judge here
                    kpoints_name_pre=[kpoints_name_pre;...
                        "K"];
                end
            end
        end
        [N_kpoints,~] = size(kpoints);
        kpoints_name_t(1)=kpoints_name_pre(1);
        kpoints_name_t(1) = strrep(kpoints_name_t(1),'GAMMA','??');
        kpoints_name_t(1) = strrep(kpoints_name_t(1),'Gamma','??');
        kpoints_name_t(1) = strrep(kpoints_name_t(1),'G','G');
        if N_kpoints >2
            for i=2:N_kpoints/2
                if strcmp(kpoints_name_pre(2*i-2,:),kpoints_name_pre(2*i-1,:))
                    kpoints_name_t(i)=kpoints_name_pre{2*i-1};
                else
                    kpoints_name_t(i)=kpoints_name_pre{2*i-2}+"|"+kpoints_name_pre{2*i-1};
                end
                kpoints_name_t(i) = strrep(kpoints_name_t(i),'GAMMA','??');
                kpoints_name_t(i) = strrep(kpoints_name_t(i),'Gamma','??');
                kpoints_name_t(i) = strrep(kpoints_name_t(i),'G','??');
            end
            kpoints_name_t(i+1)=kpoints_name_pre{2*i};
        else
            kpoints_name_t(2)=kpoints_name_pre{2};
        end
        
        % if KPOINTS_information{1} has high-symmetry information
        kpoints_string=deblank(KPOINTS_information{1});
        kpoints_string=strrep(kpoints_string,'Gamma','??');
        kpoints_string=strrep(kpoints_string,'GAMMA','??');
         kpoints_string=strrep(kpoints_string,'G','??');       
        kpoints_string=strrep(kpoints_string,' ','|');
        kpoints_name=regexp(kpoints_string, '-|', 'split');
        %if KPOINTS_information{1} has not  high-symmetry information
        if length(char(kpoints_name(1))) ~= 1 ||  length(char(kpoints_name(1))) ~= 5
            kpoints_name=kpoints_name_t;
        end
         kpoints_name(end) = strrep(kpoints_name(end),'GAMMA','??');
         kpoints_name(end) = strrep(kpoints_name(end),'Gamma','??');
         kpoints_name(end) = strrep(kpoints_name(end),'G','??');
    else
        error(fprint('Unsupport mode: \n',mode));
    end
%-------- return --------
end