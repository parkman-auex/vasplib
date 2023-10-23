%% PROCAR_read 

%% mode site

% f-mode soc-mode 
%SOC_flag=1;
%F_flag =1;
global SOC_flag F_flag ;
% bash for file operation
tic;
!head -n 3        PROCAR  > PROCAR_information
% to be extended
!grep 'k-point  ' PROCAR  > PROCAR_kpoints 
% to be extended
!grep 'band '     PROCAR  > PROCAR_eigenval
%
!sed '1,3d' PROCAR | sed '/^[ ]*$/d'|sed '/tot/d'|sed '/k-point/d'| sed '/band/d'> PROCAR_reshape 
% 3 seconds runing for a 100MB file
toc;
%% POSCAR_read
POSCAR_read;
%% read base information




% Initialize variables.
filename = 'PROCAR_information';
delimiter = ' ';
startRow = 2;
formatSpec = '%*s%*s%*s%f%*s%*s%*s%f%*s%*s%*s%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);
PROCARinformation = [dataArray{1:end-1}];
%kpoints
kpointsnum=PROCARinformation(1);
%bands
bandsnum=PROCARinformation(2);
%ions
ionsnum=PROCARinformation(3);

if ionsnum ~= sites_num
    error("wrong!!!!!,ions in POSCAR is not same with PROCAR");
end

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% PROCAR_collection
PROCAR_element=struct('ion',[],'orbital',[],'WEIGHTCAR',[],'MXWEIGHT',[],'MYWEIGHT',[],'MZWEIGHT',[]);

%ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2  fy3x2   fxyz   fyz2    fz3   fxz2   fzx2    fx3    tot
if F_flag == 1
    mode = 'f';
    orbitalsnum=17;
    PROCAR_collection=repmat(PROCAR_element,ionsnum,orbitalsnum);
%% PROCAR_reshape read 
% for a datalike(num) file matlab reads it faster a lot
    tic;
    filename = 'PROCAR_reshape';
    delimiter = ' ';
%   18 lines [f-mode]
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    PROCARreshape = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans;
    % 2 seconds for a 100MB file
    toc;


%% data split
tic       ;
% take care in SOC mode
    if SOC_flag == 1
        % the num  we use here : kpointsnum bandsnum ionsnum orbitals 
        for ion=1:ionsnum
            for orbital=1:1:orbitalsnum
                WEIGHTCAR =[];
                WEIGHTCAR_MX=[];
                WEIGHTCAR_MY=[];
                WEIGHTCAR_MZ=[];
                % WEIGHTCAR datatype
                % ************************************************************************************************************************************************
                %     ----- band     kpoint1     kpoint2     kpoint3     ............
                %           band1                                .  
                %           band2                                . 
                %             .       .     .      .      .   weight      .     .     .     .
                %             .                                  .
                %             .                                  .
                % ************************************************************************************************************************************************

                for i =1:1:bandsnum
                    for j=1:1:kpointsnum
                        % need a tool function here
                        
                        WEIGHTCAR(i,j) =PROCARreshape(PROCAR_data_location(ion,j,i,0,mode,ionsnum,bandsnum),orbital+1);
                        WEIGHTCAR_MX(i,j)=PROCARreshape(PROCAR_data_location(ion,j,i,1,mode,ionsnum,bandsnum),orbital+1);
                        WEIGHTCAR_MY(i,j)=PROCARreshape(PROCAR_data_location(ion,j,i,2,mode,ionsnum,bandsnum),orbital+1);
                        WEIGHTCAR_MZ(i,j)=PROCARreshape(PROCAR_data_location(ion,j,i,3,mode,ionsnum,bandsnum),orbital+1);
                    end
                end
        %     PROCAR_collection_data_type
        % ************************************************************************************************************************************************
        %     ----- ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2  fy3x2   fxyz   fyz2    fz3   fxz2   fzx2    fx3    tot
        %           ion1                                .  
        %           ion2                                . 
        %             .       .     .      .       WEIGHTCAR      .       .      .     .
        %             .                                 .
        %             .                                 .
        % ************************************************************************************************************************************************
                %PROCAR_element=struct('ion',[],'orbital',[],'WEIGHTCAR',[],'WEIGHTCAR_MX',[],'WEIGHTCAR_MY',[],'WEIGHTCAR_MZ',[]);        
                PROCAR_collection(ion,orbital).ion=ion;
                PROCAR_collection(ion,orbital).orbital=orbital;
                PROCAR_collection(ion,orbital).WEIGHTCAR=WEIGHTCAR;
                PROCAR_collection(ion,orbital).WEIGHTCAR_MX=WEIGHTCAR_MX;
                PROCAR_collection(ion,orbital).WEIGHTCAR_MY=WEIGHTCAR_MY;
                PROCAR_collection(ion,orbital).WEIGHTCAR_MZ=WEIGHTCAR_MZ;

            end
        end
    else
        for ion=1:ionsnum
            for orbital=1:1:orbitalsnum
                  WEIGHTCAR =[];
                % WEIGHTCAR datatype
                % ************************************************************************************************************************************************
                %     ----- band     kpoint1     kpoint2     kpoint3     ............
                %           band1                                .  
                %           band2                                . 
                %             .       .     .      .      .   weight      .     .     .     .
                %             .                                  .
                %             .                                  .
                % ************************************************************************************************************************************************

                for i =1:1:bandsnum
                    for j=1:1:kpointsnum
                        % need a tool function here                       
                        WEIGHTCAR(i,j) =PROCARreshape(PROCAR_data_location(ion,j,i,0,mode,ionsnum,bandsnum,SOC_flag),orbital+1);
                    end
                end
        %     PROCAR_collection_data_type
        % ************************************************************************************************************************************************
        %     ----- ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2  fy3x2   fxyz   fyz2    fz3   fxz2   fzx2    fx3    tot
        %           ion1                                .  
        %           ion2                                . 
        %             .       .     .      .       WEIGHTCAR      .       .      .     .
        %             .                                 .
        %             .                                 .
        % ************************************************************************************************************************************************
                %PROCAR_element=struct('ion',[],'orbital',[],'WEIGHTCAR',[],'WEIGHTCAR_MX',[],'WEIGHTCAR_MY',[],'WEIGHTCAR_MZ',[]);        
                PROCAR_collection(ion,orbital).ion=ion;
                PROCAR_collection(ion,orbital).orbital=orbital;
                PROCAR_collection(ion,orbital).WEIGHTCAR=WEIGHTCAR;
            end
        end    
    end
    %%!!!!!!!
    clearvars WEIGHTCAR WEIGHTCAR_MX WEIGHTCAR_MY WEIGHTCAR_MZ;
toc      ;
end
         
        
       
    





