%% WEIGHTCAR gen for VASPKIT -> dat
%
%
% * Label: data
%
%% Description of the Function:
%%
%% Usage: 
%
% * [WEIGHTCAR,EIGENCAR,KPATH] = WEIGHTCAR_gen(filename,seq_list,mode) 
% * [WEIGHTCAR,EIGENCAR,KPATH] = WEIGHTCAR_gen(filename,seq_list) 
% * [WEIGHTCAR,EIGENCAR,KPATH] = WEIGHTCAR_gen(filename) 
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
% * Document Date: 2020/12/04
% * Creation Date: 2020/12/04
% * Last updated : 2020/12/04
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code : 
%
%
function [WEIGHTCAR,EIGENCAR,KPATH] = WEIGHTCAR_gen(filename,seq_list,mode) 
%--------  init  --------
   import vasplib_tool.*
   import linux_matlab.*
%--------  narg  --------
    if nargin < 3
        mode = 'vaspkit';
    end
    if nargin < 2
        seq_list = -1 ;
    end
%--------  chek  --------

%--------  fbug  --------

%--------  juge  --------
    if strcmp(mode,'vaspkit-band')
        fprintf('support vaspkit 1.2.0');
        %         command_str = "!cp "+filename+" temp_Pband.dat";
        %         eval(command_str);
        % copy for delete
        %copyfile(filename,'temp_Pband.dat');
        [Bandindex,rm_num_list] = grep(filename,'#','silence');

        rm_line(filename,rm_num_list,'temp_Pband.dat');
        WEIGHTCAR_init = load('temp_Pband.dat');
        %--------  init  --------
        Band_index = length(Bandindex)-2;
        [TOTnum,width] = size(WEIGHTCAR_init);
        K_num = TOTnum/Band_index;
        %--------  print  --------
        if width >10
            fprintf('s py pz px dxy dyz dz2 dxz dx2-y2 fy3x2  fxyz  fyz2 fz3 fxz2  fzx2  fx3 tot');
            fprintf('1  2  3  4   5   6   7   8      9    10    11    12  13   14  	 15   16  17');
        else
            fprintf('s py pz px dxy dyz dz2 dxz dx2-y2 tot');
            fprintf('1  2  3  4   5   6   7   8      9  10');
        end
        KPATH = reshape(WEIGHTCAR_init(:,1),K_num,Band_index);
        KPATH = KPATH(:,1);
        EIGENCAR = reshape(WEIGHTCAR_init(:,2),K_num,Band_index);
        EIGENCAR = EIGENCAR.';
        WEIGHTCAR_init(:,[1,2]) = [];
        if seq_list > 0
            WEIGHTCAR_init2 = sum(WEIGHTCAR_init(:,seq_list),2);
            WEIGHTCAR = reshape(WEIGHTCAR_init2,K_num,Band_index).';
        else
            disp('WEIGHTCAR with whole projection');
            for i =1:width-2
                WEIGHTCAR(:,:,i) = reshape(WEIGHTCAR_init(:,i),K_num,Band_index).';
            end
        end
    elseif strcmp(mode,'vaspkit-band-silence')
        [Bandindex,rm_num_list] = grep(filename,'#','silence');        
        rm_line(filename,rm_num_list,'temp_Pband.DAT');
        WEIGHTCAR_init = load('temp_Pband.DAT');
        %--------  init  --------
        Band_index = length(Bandindex)-2;
        [TOTnum,width] = size(WEIGHTCAR_init);
        K_num = TOTnum/Band_index;
        %--------  print  --------
        KPATH = reshape(WEIGHTCAR_init(:,1),K_num,Band_index);
        KPATH = KPATH(:,1);
        EIGENCAR = reshape(WEIGHTCAR_init(:,2),K_num,Band_index);
        EIGENCAR = EIGENCAR.';
        WEIGHTCAR_init(:,[1,2]) = [];
        if seq_list > 0
            WEIGHTCAR_init2 = sum(WEIGHTCAR_init(:,seq_list),2);
            WEIGHTCAR = reshape(WEIGHTCAR_init2,K_num,Band_index).';
        else
            for i =1:width-2
                WEIGHTCAR(:,:,i) = reshape(WEIGHTCAR_init(:,i),K_num,Band_index).';
            end
        end
    elseif strcmp(mode,'vaspkit-DOS')
        disp('support vaspkit 1.2.0 - DOS');
        WEIGHTCAR_init = textread(filename,'','headerlines',1);
        [~,width] = size(WEIGHTCAR_init);
        if width < 12
            disp('  Energy     s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot') ;
            disp('s py pz px dxy dyz dz2 dxz dx2-y2 tot');
            disp('1 2  3  4  5   6   7   8   9  	10');
        elseif width >12
            disp('   Energy 	    s  	   py  	   pz  	   px  	  dxy  	  dyz  	  dz2  	  dxz  	x2-y2  	fy3x2  	 fxyz  	 fyz2  	  fz3  	 fxz2  	 fzx2  	  fx3  	  tot');
            disp('s py pz px dxy dyz dz2 dxz dx2-y2 fy3x2  fxyz  fyz2 fz3 fxz2  fzx2  fx3 tot');
            disp('1 2  3  4  5   6   7   8   9  	10     11  	 12   13  14  	15    16  17');
        end
        EIGENCAR = WEIGHTCAR_init(:,1);
        WEIGHTCAR_init(:,1) = [];
        if seq_list > 0
            WEIGHTCAR = sum(WEIGHTCAR_init(:,seq_list),2);
        else
            disp('WEIGHTCAR with whole projection');
            for i =1:width-2
                WEIGHTCAR = WEIGHTCAR_init;
            end
        end
    elseif strcmp(mode,'vaspkit-DOS-silence')
        disp('support vaspkit 1.2.0 - DOS');
        KPATH = [];
        WEIGHTCAR_init = textread(filename,'','headerlines',1);
        [TOTnum,width] = size(WEIGHTCAR_init);
        EIGENCAR = WEIGHTCAR_init(:,1);
        WEIGHTCAR_init(:,1) = [];
        if seq_list > 0
            WEIGHTCAR = sum(WEIGHTCAR_init(:,seq_list),2);
        else
            WEIGHTCAR = WEIGHTCAR_init;
        end
    elseif strcmp(mode,'PROCAR')
    elseif strcmp(mode,'WT')
    elseif strcmp(mode,'origin file')
    end

%-------- return --------    
end


