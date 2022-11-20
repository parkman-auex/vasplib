%% To Get EIGENCAR from file EIGENVAL
%
%
%% Description of the Function:
%
%% Usage: 
%
% * [EIGENCAR,EIGENCAR2,Efermi]=EIGENVAL_read(mode,EIGENVAL,Efermi)
% * [EIGENCAR,EIGENCAR2,Efermi]=EIGENVAL_read(mode,EIGENVAL)
% *  [EIGENCAR,EIGENCAR2,Efermi]=EIGENVAL_read(mode)
%
%% Input:
%  
% # Efermi: for handy
% # EIGENVAL: filename
% # mode: vasp qe vaspkit
%
%% Output:
%
% # EIGENCAR :
% # EIGENCAR2: for spinful
% # Efermi :
%
%% example:
%%% for VASP_spinless/SOC:
%   EIGENCAR=EIGENVAL_read()
%%% for VASP_spinful: 
%   [EIGENCAR,EIGENCAR] = EIGENVAL_read()
%%% for vaspkit:
%   EIGENCAR = EIGENVAL_read('vaspkit','BAND.dat',0)
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
function [EIGENCAR,EIGENCAR2,Efermi]=EIGENVAL_read(mode,EIGENVAL,Efermi)
    %--------  init  --------
    %I = 1i;
    %--------  narg  --------
    if nargin <1
        mode = 'vasp';
    end
    if nargin < 2
        if strcmp(mode,'vasp')
            EIGENVAL = 'EIGENVAL';
        elseif strcmp(mode,'qe')
            EIGENVAL = 'BAND.dat';
        else
            EIGENVAL = 'BAND.dat';
        end
    end
    if nargin <3
        if strcmp(mode,'vaspkit')
            Efermi = 0;
        elseif exist('Efermi','file')
            Efermi = double(textscan('Efermi'));
        elseif exist('DOSCAR','file') && strcmp(mode,'vasp')
            Efermi = GetFermi('vasp');
        elseif exist('scf.out','file') && strcmp(mode,'qe')
            Efermi = GetFermi('qe');
        else
            Efermi = 0;
        end
    end
    %--------  chek  --------
    if strcmp(mode,'vasp')
        % read EIGENVAL
        data=textread(EIGENVAL,'','headerlines',5);
        NBands=data(1,3);
        Nelectrons=data(1,1);
        Ktotal=data(1,2);
        % setup EIGENCAR
        data=textread(EIGENVAL,'','headerlines',8);
        for i=1:1:NBands
            for j=1:1:Ktotal
                Pup(j,: )=[j,data((NBands+1)*(j-1)+i,2)-Efermi];         %
                Pdown(j,: )=[j,data((NBands+1)*(j-1)+i,3)-Efermi];         %
            end
            EIGENCAR(i,:)=Pup(:,2)';
            EIGENCAR2(i,:)=Pdown(:,2)';
        end
    elseif strcmp(mode,'qe')
        % import linux_matlab.*;
        [~,num_list] = linux_matlab.grep(EIGENVAL,'    0.0000','silence');
        NBands = length(num_list)  ;
        data = load(EIGENVAL)    ;
        [tot_row,~] = size(data );
        Ktotal = tot_row/NBands  ;
        EIGENCAR = reshape(data(:,2),Ktotal,NBands);
        EIGENCAR = EIGENCAR.' - Efermi*ones(NBands,Ktotal);
    elseif strcmp(mode,'vaspkit')
        % import linux_matlab.*;
        [~,num_list] = linux_matlab.grep(EIGENVAL,'    0.0000','silence');
        NBands = length(num_list)  ;
        data = importdata(EIGENVAL)    ;
        [tot_row,~] = size(data );
        Ktotal = tot_row/NBands  ;
        EIGENCAR = reshape(data(:,2),Ktotal,NBands);
        EIGENCAR = EIGENCAR.' - Efermi*ones(NBands,Ktotal);
    end

end
