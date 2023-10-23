function [ PARCHGCAR , PARCHG_line ]= PARCHG_gen(orb_list,WaveFunc,PARCHG_file,options)
% PARCHG gen
% usage: PARCHG_gen(orb_list,WaveFunc,PARCHG_file,curve_R,Rm,sites,Atom_name,Atom_num )
%        PARCHG_gen(orb_list,WaveFunc,PARCHG_file,curve_R)
%    a    PARCHG_gen(orb_list,WaveFunc,PARCHG_file)
%        PARCHG_gen(orb_list,WaveFunc)
% PARCHG contains the partial charge densities.
% This file contains
%   the lattice vectors,
%	atomic coordinates,
%	the total charge density multiplied by the volume { \rho (r)*V_{{{\rm {cell}}}}}\rho (r)*V_{{{\rm {cell}}}} on the fine FFT-grid (NG(X,Y,Z)F),
% and the PAW one-center occupancies.
% WRITE(IU,FORM) (((C(NX,NY,NZ),NX=1,NGXF),NY=1,NGYF),NZ=1,NGZF)
% Note that the real-space mesh (NX,NY,NZ) is uniform and is spanned by the primitive lattice vectors
% (a,b,c)defined in the POSCAR and can read explicitly
% \left(N_{x}, N_{y}, N_{z}\right) \triangleq \frac{N_{x}-1}{N_{G X P}} \vec{a}+\frac{N_{y}-1}{N_{G Y F}} \vec{b}+\frac{N_{z}-1}{N_{G Z F}} \vec{c}
%% naigin
arguments
    orb_list ;
    WaveFunc ;
    PARCHG_file = 'PARCHG.vasp';
    options.cutoff = 1.53 ;
    options.unit = [0.2,0.2,0.2]; % divide your lattice
    options.Rm = [];
    options.sites = [];
    options.Atom_name = [];
    options.Atom_num  = [];
    options.ax =  handle([]);
    options.POSCAR = 'POSCAR';
    options.WaveMin = 1e-3;
    options.WaveColor = 'r';
    options.WaveSize = 1;
    options.OrbColor = 'k';
    options.OrbSize = 1;
end
if nargin < 3
    mode = 'glance';
else
    mode = 'PARCHG';
end
if nargin <2
    WaveFunc = zeros(size(orb_list,1),1);
end
if isempty(options.sites) && isempty(options.Atom_name) && isempty(options.Atom_num)
    [Rm,sites,Atom_name,Atom_num,~]=POSCAR_read(options.POSCAR );
else
    Rm = options.Rm;
    sites = options.sites;
    Atom_name = options.Atom_name;
    Atom_num = options.Atom_num;
end
[Nwave,Nlist] = size(WaveFunc);
if strcmp(mode,'glance')
    waveplot(orb_list,WaveFunc);
elseif strcmp(mode,'PARCHG')
    %% check
    Norb = length(orb_list);
    disp(Norb);
    if  Norb ~= Nwave
        error('Orbital list length is not equal to WaveFunc');
    end
    %% init
    %     adjustparm = 30;
    if nargin <4
        curve_R = [0.2 0.2 0.2];
    end
    %adjustparm = round(3/norm(curve_R))+1;
    
    a = norm(Rm(1,:));
    b = norm(Rm(2,:));
    c = norm(Rm(3,:));
    % Rmesh init
    NGXF = round(a/curve_R(1));
    NGYF = round(b/curve_R(2));
    NGZF = round(c/curve_R(3));
    PARCHGCAR = zeros(NGXF,NGYF,NGZF);
    adjustparmX = round(2000/NGXF);
    adjustparmY = round(2000/NGYF);
    adjustparmZ = round(2000/NGZF);
    %     disp(NGXF);
    %     disp(NGYF);
    %     disp(NGZF);
    %% Data
    % init PARCHG_line
    PARCHG_line_single = struct('seq',[],'NX',[],'NY',[],'NZ',[],'Rc',[],...
        'To_orbital_list',[],'To_orbital_seq_list',[],'To_orbital_Rlength_list',[],'To_orbital_Wave_amp_list',[],'addornot',0);
    PARCHG_line = repmat(PARCHG_line_single,[NGXF*NGYF*NGZF,1]);
    
    tic;
    for iorb = 1:Norb
        temp_orb = (orb_list(iorb,:));
        Ro(1) = round(temp_orb(1)*NGXF);
        Ro(2) = round(temp_orb(2)*NGYF);
        Ro(3) = round(temp_orb(3)*NGZF);
        %Ro = round(temp_orb./curve_R);
        %disp(Ro);
        
        for kz = Ro(3)-adjustparmZ : Ro(3)+adjustparmZ+1
            if kz >=1 && kz <= NGZF
                for jy = Ro(2)-adjustparmY : Ro(2)+adjustparmY+1
                    if jy >=1 && jy <= NGYF
                        for ix = Ro(1)-adjustparmX : Ro(1)+adjustparmX+1
                            if ix >=1 && ix <= NGXF
                                seq = (kz-1)*NGYF*NGXF + (jy-1)*NGXF + ix;
                                %disp(kz);
                                %disp(seq);
                                PARCHG_line(seq).seq = seq;
                                PARCHG_line(seq).NX  = ix;
                                PARCHG_line(seq).NY  = jy;
                                PARCHG_line(seq).NZ  = kz;
                                PARCHG_line(seq).Rc  = [(ix-1)/NGXF,(jy-1)/NGYF,(kz-1)/NGZF];
                                %disp(PARCHG_line(seq).Rc);
                                [add_orb,Rlength] = check_orbital_inrange(PARCHG_line(seq).Rc,temp_orb,Rm,options.cutoff);
                                if add_orb == 1
                                    %disp('bug_fix');
                                    PARCHG_line(seq).To_orbital_list = [PARCHG_line(seq).To_orbital_list;temp_orb];
                                    PARCHG_line(seq).To_orbital_seq_list=[PARCHG_line(seq).To_orbital_seq_list;iorb];
                                    PARCHG_line(seq).To_orbital_Rlength_list=[PARCHG_line(seq).To_orbital_Rlength_list;Rlength];
                                    temp_C = 0;
                                    for j =1:Nlist
                                        temp_amp = WaveFunc(iorb,j);
                                        temp_C=temp_C+temp_amp*temp_amp'*1000;
                                    end
                                    WFamp = temp_C;
                                    PARCHG_line(seq).To_orbital_Wave_amp_list=[PARCHG_line(seq).To_orbital_Wave_amp_list;WFamp];
                                    PARCHG_line(seq).addornot = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    toc;
    
    tic;
    cont = 0;
    for k = 1:NGZF
        for j = 1:NGYF
            for i = 1:NGXF
                cont = cont +1 ;
                seq = (k-1)*NGYF*NGXF + (j-1)*NGXF + i;
                if seq ~= cont
                    disp(cont);
                    disp(seq);
                    error('!!!!');
                end
                % struct
                if PARCHG_line(cont).addornot == 1
                    PARCHG_line(cont).C = S_wave_amp(PARCHG_line(cont).To_orbital_Rlength_list,PARCHG_line(cont).To_orbital_Wave_amp_list) ;
                    PARCHGCAR(i,j,k) =   PARCHG_line(cont).C ;
                else
                    PARCHG_line(cont).seq = cont;
                    PARCHG_line(cont).NX  = i;
                    PARCHG_line(cont).NY  = j;
                    PARCHG_line(cont).NZ  = k;
                    PARCHG_line(cont).Rc  = [(i-1)/NGXF,(j-1)/NGYF,(k-1)/NGZF];
                    PARCHG_line(cont).C = 0;
                    PARCHGCAR(i,j,k) =   PARCHG_line(cont).C ;
                end
            end
        end
    end
    toc;
    %% Write
    tic;
    %disp(PARCHG_file);
    fileID = fopen(PARCHG_file,'w');
    % init
    a_crystal_constance=1;
    title = "PARCHG Gen by MATLAB (parkman) ";
    % -------------------- write POSCAR first --------------------
    fprintf(fileID,"%s\n",title);
    fprintf(fileID,"%d\n",a_crystal_constance);
    %fprintf(fileID,"  ",Rm(i,j));
    for i=1:3
        for j=1:3
            fprintf(fileID,"  %f",Rm(i,j));
        end
        fprintf(fileID,"\n");
    end
    for i=1:length(Atom_name)
        fprintf(fileID,"%s ",Atom_name(i));
    end
    fprintf(fileID,"\n  ");
    for i=1:length(Atom_num)
        fprintf(fileID,"%d ",Atom_num(i));
    end
    fprintf(fileID,"\n");
    fprintf(fileID,"Direct\n  ");
    % sites
    [~,sites_num]=size(sites);
    for i=1:sites_num
        fprintf(fileID,"%f  ",sites(i).rc1);
        fprintf(fileID,"%f  ",sites(i).rc2);
        fprintf(fileID,"%f  ",sites(i).rc3);
        %         if ~strcmp(string(sites(i).name),"")
        %             %fprintf(fileID,"%s\n  ",sites(i).name);
        %             fprintf(fileID,"\n  ");
        %         else
        %             fprintf(fileID,"\n  ");
        %         end
        fprintf(fileID,"\n  ");
    end
    fprintf(fileID,"\n  ");
    % -------------------- write Nx Ny Nz GF second --------------------
    fprintf(fileID,"%d %d %d \n ",NGXF,NGYF,NGZF);
    % -------------------- write Nx Ny Nz GF second --------------------
    Ten_cont = 0;
    for k = 1:NGZF
        for j = 1:NGYF
            for i = 1:NGXF
                %seq = (k-1)*NGYF + (j-1)*NGXF + i;
                
                Ten_cont = Ten_cont + 1;
                %fprintf(fileID,"%8.4f  ", PARCHG_line(seq).C);
                fprintf(fileID,"%8.6f  ", PARCHGCAR(i,j,k));
                if Ten_cont == 10
                    Ten_cont = 0;
                    fprintf(fileID,"\n");
                end
            end
        end
    end
    fclose(fileID);
    toc;
end
end

function [checkvalue,lengthR12]= check_orbital_inrange(Rc,orb,Rm,cutoff)
    if nargin < 5
        cutoff  = 1.53;
    end

    % To_orbital_list = [];
    % To_orbital_seq_list = [];
    % To_orbital_Rlength_list = [];
    % To_orbital_Wave_amp_list = [];
    lengthR12 = Rlength(Rc,orb,Rm);
    %         disp(lengthR12);
    if lengthR12 < cutoff
        checkvalue = 1;
        %             To_orbital_list = [To_orbital_list;orb_list(i,:)];
        %             To_orbital_seq_list = [To_orbital_seq_list;i];
        %             To_orbital_Rlength_list =[To_orbital_Rlength_list;lengthR12];
        %             To_orbital_Wave_amp_list = [To_orbital_Wave_amp_list;WaveFunc(i)];
    else
        checkvalue = 0;
    end

end



function C = S_wave_amp(Rlength_list,Wave_amp_list)
    % nargin
    % if nargin < 3
    %     adjustparm = 20;
    % end
    % check
    [NWave_amp,Nlist] = size(Wave_amp_list);
    % if NWave_amp ~= Rlength_list
    %     error('nRlength_list ~= nWave_amp, why??');
    % end
    C_list = [];
    for i =1:NWave_amp
        Rr = Rlength_list(i,:);
        if Rr <  1

            Rr =1;
        end
        temp_C=0;
        for j =1:Nlist
            temp_amp = Wave_amp_list(i,j);
            if temp_amp*temp_amp' < 1e-4
                temp_amp =0;
            end
            temp_C=temp_C+temp_amp*temp_amp'*exp(-Rr/0.53)^2;
        end
        C_list = [C_list ;temp_C];
    end
    C = sum(C_list);

end

function lengthR12 = Rlength(Rc1,Rc2,Rm)
    lengthR12 = norm((Rc2-Rc1)*Rm);
end
