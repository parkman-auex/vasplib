%% PROCAR_data_location
% NUMFLAG :
% 0       : TOT
% 1       : MX
% 2       : MY
% 3       : MZ
function line = PROCAR_data_location(ion,kpoint,band,NUMFLAG,mode,ionsnum,bandsnum,SOC_flag)
    if nargin <8
        SOC_flag =1;
    end
    if nargin < 6
        !head -n 3        PROCAR  > PROCAR_information
        % Initialize variables.
        filename = 'PROCAR_information';
        delimiter = ' ';
        startRow = 2;
        formatSpec = '%*s%*s%*s%f%*s%*s%*s%f%*s%*s%*s%f%[^\n\r]';
        fileID = fopen(filename,'r');
        textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false);
        fclose(fileID);
        PROCARinformatiobn = [dataArray{1:end-1}];
        %kpoints
        kpointsnum=PROCARinformatiobn(1);
        %bands
        bandsnum=PROCARinformatiobn(2);
        %ions
        ionsnum=PROCARinformatiobn(3);
        % Clear temporary variables
        clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    end
    if nargin < 5
        mode = 'f' ;
    end
     % to be extended
%%
    if strcmp(mode,'f')
        % to be extended
        if SOC_flag == 1;
            kpoints_round = bandsnum*ionsnum*4;
            bands_round   = ionsnum*4;
            line = (kpoint-1)*kpoints_round + (band-1)*bands_round + NUMFLAG*ionsnum+ion;
        else
            kpoints_round = bandsnum*ionsnum;
            bands_round   = ionsnum;
            line = (kpoint-1)*kpoints_round + (band-1)*bands_round + NUMFLAG*ionsnum+ion;
        end
    end

end