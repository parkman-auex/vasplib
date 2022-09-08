%% Read wannier90_hr.dat
% Be careful of the struct
% usage [H_xyz,row] =Readhr(hrdat)
% usage [H_xyz,row] =Readhr(hrdat,mode) mode = 'n'
% H_xyz_ =struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);

function [H_xyz,row] = Read_hr(hrdat,mode)
    tic;
    if nargin <2
        mode = 'n';
    elseif size(mode) == [3,3]
        mode = 'n';
    else
        mode = 'sym';
        syms Rm a11 a12 a13 a21 a22 a23 a31 a32 a33 real;
        Rm = [a11 a12 a13;a21 a22 a23;a31 a32 a33];
        disp("Model Mode");
    end
    
    %% read hopping terms
    if nargin < 1
        hrdat = 'wannier90_hr.dat';
    end

    %% read Nbands and NRPTS
    delimiter = ' ';
    startRow = 2;
    endRow = 3;
    formatSpec = '%d%*s%[^\n\r]';
    fileID = fopen(hrdat,'r');
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter,...
        'MultipleDelimsAsOne', true, 'TextType', 'string',...
        'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    NBR = dataArray{:, 1};
    NUM_WAN = NBR(1);
    NRPTS = NBR(2);
    NRPTS_num1=fix(double(NRPTS)/15);
    NRPTS_num2=mod(NRPTS,15);
    fclose(fileID);
    %% read NRTT_list
    % NRPTS_num1=floor(NRPTS,15);
    % NRPTS_num2=rem(NRPTS,15);
    fileID = fopen(hrdat,'r');
    startRow = 4;
    endRow = startRow+NRPTS_num1;% caculate here
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec,NRPTS_num1+1 , 'Delimiter', delimiter,...
        'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', 0,...
        'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    NRPT_list = [dataArray{1:end-1}];
    NRPT_list = reshape(NRPT_list',NRPTS_num1*15+15,1);
    fclose(fileID);
    %% read hopping terms
    fileID = fopen(hrdat,'r');
    if NRPTS_num2==0
        startRow = endRow;
    else
        startRow = endRow+1;
    end

%     formatSpec = '%f%f%f%f%f%f%f%*s%[^\n\r]';
%     dataArray = textscan(fileID, formatSpec, 'Delimiter', ' ', 'WhiteSpace', '',...
%         'HeaderLines' ,startRow-1, 'ReturnOnError', false);
%     fclose(fileID);
    formatSpec = '%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';


    %% 根据格式读取数据列。
    % 该调用基于生成此代码所用的文件的结构。如果其他文件出现错误，请尝试通过导入工具重新生成代码。
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    %% 关闭文本文件。
    fclose(fileID);

    Vec_Fir = dataArray{:, 1};      % a1 direction
    Vec_Sec = dataArray{:, 2};      % a2 direction
    Vec_Thi = dataArray{:, 3};      % a3 direction
    Orb_fir = dataArray{:, 4};
    Orb_sec = dataArray{:, 5};
    h_real = dataArray{:, 6};
    h_imag = dataArray{:, 7};
    %disp(h_imag);
    clearvars filename startRow formatSpec fileID dataArray ans;
    %% Clear temporary variables
    % read

    V_f = reshape(Vec_Fir, NUM_WAN*NUM_WAN, NRPTS);
    V_s = reshape(Vec_Sec, NUM_WAN*NUM_WAN, NRPTS);
    V_t = reshape(Vec_Thi, NUM_WAN*NUM_WAN, NRPTS);
    V=[V_f(1,:)',V_s(1,:)',V_t(1,:)'];
    [~,row]=ismember([0 0 0],V,'rows');
    h_r = reshape(h_real, NUM_WAN*NUM_WAN, NRPTS);
    h_i = reshape(h_imag, NUM_WAN*NUM_WAN, NRPTS);

    H_xyz_ =struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
    H_xyz = repmat(H_xyz_ ,[NRPTS,1]);    % Hamiltonian of every cell;

    syms k_x k_y k_z real
    if strcmp(mode,'n')
        for n = 1 : NRPTS
            H_xyz(n).seq=n;
            H_xyz(n).Degen=NRPT_list(n);
            H_xyz(n).Hnum =  reshape(h_r(:,n), NUM_WAN, NUM_WAN) + 1i*reshape(h_i(:,n), NUM_WAN, NUM_WAN);
            H_xyz(n).Hcoe = H_xyz(n).Hnum;
            H_xyz(n).vector = V(n,1:3);
            %From here , get hints
        end
    else
        %% you need a function here
        %H_xyz(n).key
        %H_xyz(n).nokey
        %H_xyz(n).Hsym
        %H_xyz(n).Hstr
        H_xyz(n).Hnum =  reshape(h_r(:,n), NUM_WAN, NUM_WAN) + 1i*reshape(h_i(:,n), NUM_WAN, NUM_WAN);
        H_xyz(n).Hcoe = sym(H_xyz(n).Hnum);
        H_xyz(n).vector = V(n,1:3);
        H_xyz(n).Hsym = H_xyz(n).Hcoe * exp(1i*[k_x k_y k_z]*Rm'*sym(H_xyz(n).vector)');
        H_xyz(n).Hfun=matlabFunction(H_xyz(n).Degen*H_xyz(n).Hsym,'Vars',{k_x,k_y,k_z});
    end
    %    H_xyz_0 = H_xyz(row);
    %    H_xyz(row)=[];
    %    H_xyz=[H_xyz_0;H_xyz];





    toc;
end









