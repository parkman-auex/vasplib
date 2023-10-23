%Read HSPICE generated ASCII formatted .tr#/.sw#/.ac# files 
%% Author
%Mohammad Abu Raihan Miah
%University of California, San Diego
%ver 4.0.0, 11/23/21
%% Function description
%This function reads a ASCII formatted HSPICE output file (.option post=2) 
%'filename.tr#' or 'filename.sw#' or 'filename.ac#' and saves all the 
%signals in the variable 'simulation_result' as a structure. 
%
%simulation_result=a structure (2 fields) that contains the contents from
%                  the HSPICE file.  
%simulation_result(#).var_name=name of the signals present in the file.
%simulation_result(#).val=a vector contains the values of the signal with
%                         the name simulation_result(#).var_name
%
%% Example for calling this function:
% sim_data=read_hspice_tr_sw_ac('tr_example.tr0');
% sim_data=read_hspice_tr_sw_ac('sw_example.sw0');
% sim_data=read_hspice_tr_sw_ac('ac_example.ac0');
%% main function
%Don't touch here
function simulation_result=read_hspice_tr_sw_ac(filename)

    % determine the file extension
    [filepath,name,ext] = fileparts(filename);clear filepath name
    file_extension=ext(1:3); clear ext

    % read data from the file
    file_handle = fopen(filename);
    content_file = textscan(file_handle,'%s','HeaderLines',3);
    fclose(file_handle);

    % find number of the variables and start point of the data
    content_str=string(content_file{1,1});clear content_file
    data_start_ind=find(content_str==string('$&%#'))+1;
    num_var=find(isnan(str2double(content_str(1:data_start_ind-2))),1)-1;
    
    % save name of the signals in the .var_name field
    save_signal_names();
    clear num_var

    % save signal data in the .val field
    data=replace(strjoin(content_str(data_start_ind:end)),' ','');
    data_frmt=char(ones(1,find(char(content_str(data_start_ind))=='E',1)-1)*'.');
    data_separate=double(regexp(data,[data_frmt 'E[+-]..'],'match'));clear data
    for ii=1:length(simulation_result)
        simulation_result(ii).val=data_separate(ii:length(simulation_result):length(data_separate)-1)';
    end

    function save_signal_names()
        var_name_raw_ind=num_var+2:data_start_ind-2;
        temp=find(max(max((char(content_str(var_name_raw_ind))=='(')*2,1),[],2)>1);
        var_name_1st_part=content_str(num_var+1+temp);
        var_name_special=setdiff(var_name_raw_ind,var_name_raw_ind(temp));
        var_name_last_part(1:length(temp),1)=string('');
        for iii=1:length(var_name_special)
            tt=find(~(var_name_special(iii)>var_name_raw_ind(temp)),1)-1;
            var_name_last_part(tt)=content_str(var_name_raw_ind(tt+1));
        end
        var_name_raw=[content_str(num_var+1) var_name_1st_part'+var_name_last_part'+string(')')]';
        
        % add _real and _imag with the relevant variable names from .ac#
        % files
        if file_extension == '.tr' | file_extension == '.sw'
            simulation_result=cell2struct(cellstr(var_name_raw'),'var_name',1);
        elseif file_extension == '.ac'
            var_real_imag=find(rem(str2double(content_str(1:num_var)),7)==1);
            if isempty(var_real_imag)
                var_name=var_name_raw;
            else
                count=1;
                for iii=1:length(var_name_raw)
                    if any(iii==var_real_imag)
                        var_name(count)=var_name_raw(iii)+string('_real');count=count+1;
                        var_name(count)=var_name_raw(iii)+string('_imag');
                    else
                        var_name(count)=var_name_raw(iii);
                    end
                    count=count+1;
                end
            end
            simulation_result=cell2struct(cellstr(var_name),'var_name',1);
        else
            error(['The toolbox cannot handle' file_extension '# files']);
        end
    end
end
