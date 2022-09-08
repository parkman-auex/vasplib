%% Gen wannier90_hr.dat
% Be careful of the struct

% parkman 2019.9.23 fisrt
% parkman 2020.3.07 second

%% properties

% input :
%       H_xyz      (To Hr.dat data format)
        %H_xyz_ =struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]); 
%       wann_num (dimention of matrix)
%       NRPTS 
%       NRPTS_list[]      (NRPT_degeneracy sequence) 
%       \\R_abc    (3 dimension lattice matrix)
% output :
%       hr.dat

%% format of hr.dat

%1 date and time 
%2 [tab]wann_num (same with .win file)
%3 [tab]nprt     ( in .win file :  )
%4 degeneracy of each Grid;15 num per line!!!!! 
%(15mod$nprt+1+5) Rk(1 2 3) M_m M_n H_r(real_part) H_i(imag_part)

%% function may be used
% 
% each H attain M_H and Traverse every elements get its real and imaginary part
% formational output

%% main_function
function [hrdat,splen]= Gen_hr(H_xyz,output,mode,infinit_small)
tic;
splen =0;
%% read hopping terms
if nargin < 2
    ouput = 'wannier90_hr.dat';
end
if nargin < 3
    mode = 'dense';
end
if nargin < 4
    infinit_small = 1e-8;
end

hrdat=fopen(output,'w');
[NRPTS,~]=size(H_xyz);
WAN_NUM=length(H_xyz(1).Hnum);

if strcmp(mode,'dense')
    %% write title Nbands and NRPTS
        date_=date;
        fprintf(hrdat,"%s\n",date_);
    %1
    %2
        fprintf(hrdat,"         %d\n",WAN_NUM);
    %2
    %3
        fprintf(hrdat,"         %d\n",NRPTS);
    %3
    %4
    %% write NRPT_list
    % NRPTS_num1=fix(NRPTS/15);
    % NRPTS_num2=rem(NRPTS,15);
    COUNT_FLAG=0;
    for i =1:NRPTS
        fprintf(hrdat,"    %d",H_xyz(i).Degen);
        COUNT_FLAG=COUNT_FLAG+1;
        if COUNT_FLAG == 15
            COUNT_FLAG = 0;
            fprintf(hrdat,"\n");
        end
    end
    %another \n
    fprintf(hrdat,"\n");
    %% write hopping terms
        for i=1:NRPTS
            fprintf('Wrinting £¨%d/%d £© NRPT\n',i,NRPTS);
            for k=1:WAN_NUM
                %fprintf('Wrinting %d th WAN_ORB \n',k);
                for j=1:WAN_NUM           
                        fprintf(hrdat,"    %d    %d    %d",H_xyz(i).vector(1),H_xyz(i).vector(2),H_xyz(i).vector(3));
                        fprintf(hrdat,"    %d    %d",j,k);
                        real_part=real(H_xyz(i).Hnum(j,k));
                        imag_part=imag(H_xyz(i).Hnum(j,k));
                        fprintf(hrdat,"    %.7f    %.7f\n",real_part,imag_part);               
                end
            end
        end
    toc;
elseif strcmp(mode,'sparse')
    %% write title Nbands and NRPTS
        date_=date;
        fprintf(hrdat,"%s\n",date_);
    %1
        fprintf(hrdat,"%s\n",'splen');
        fprintf(hrdat,"         %d\n",WAN_NUM);
    %2
    %3
        fprintf(hrdat,"         %d\n",NRPTS);
    %3
    %4
    %% write NRPT_list
    % NRPTS_num1=fix(NRPTS/15);
    % NRPTS_num2=rem(NRPTS,15);
    COUNT_FLAG=0;
    for i =1:NRPTS
        fprintf(hrdat,"    %d",H_xyz(i).Degen);
        COUNT_FLAG=COUNT_FLAG+1;
        if COUNT_FLAG == 15
            COUNT_FLAG = 0;
            fprintf(hrdat,"\n");
        end
    end
    %another \n
    fprintf(hrdat,"\n");
    splen =0;
    %% write hopping terms
        for i=1:NRPTS
            fprintf('Wrinting £¨%d/%d £© NRPT\n',i,NRPTS);
            for k=1:WAN_NUM
                %fprintf('Wrinting %d th WAN_ORB \n',k);
                for j=1:WAN_NUM    
                    if abs(norm(H_xyz(i).Hnum(j,k))) >  infinit_small
                        splen =splen+1;
                        fprintf(hrdat,"    %d    %d    %d",H_xyz(i).vector(1),H_xyz(i).vector(2),H_xyz(i).vector(3));
                        fprintf(hrdat,"    %d    %d",j,k);
                        real_part=real(H_xyz(i).Hnum(j,k));
                        imag_part=imag(H_xyz(i).Hnum(j,k));
                        fprintf(hrdat,"    %.7f    %.7f\n",real_part,imag_part);   
                    end
                end
            end
        end
    toc;     
end

end


