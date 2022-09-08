function [ PARCHGCAR  ]= PARCHG_gen_wire(orb_list,WaveFunc,Ncycle,PARCHG_file,curve_R,Rm,sites,Atom_name,Atom_num )
% PARCHG gen
% usage: PARCHG_gen_wire(orb_list,WaveFunc,PARCHG_file,curve_R,Rm,sites,Atom_name,Atom_num )
%        PARCHG_gen_wire(orb_list,WaveFunc,PARCHG_file,curve_R)
%        PARCHG_gen_wire(orb_list,WaveFunc,PARCHG_file)
%        PARCHG_gen_wire(orb_list,WaveFunc)
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
if nargin < 3 
    Ncycle = 1;
end
if nargin < 4
    disp('glance mode')
    %PARCHG_file = "PARCHG";
    mode = 'glance';
else
    mode = 'PARCHG';
end
if nargin < 5
   curve_R = [0.2 0.2 0.2]; 
end

if nargin < 6
    if exist('POSCAR','file')
        [Rm,sites,Atom_name,Atom_num,~]=POSCAR_readin('POSCAR','vasp');
    else
        error('POSCAR need！');
    end
end
[Nwave,Nlist] = size(WaveFunc);

% HSVCAR = HSVCAR_gen(orb_list,'hinge');
% WEIGHT_affix =  abs(-HSVCAR(:,1))*4;
% WaveFunc_new = WaveFunc.*WEIGHT_affix;


if strcmp(mode,'glance')
%     HSVCAR = HSVCAR_gen(orb_list,'hinge');
%     WEIGHT_affix =  abs(-HSVCAR(:,1))*4;
    % WaveFunc = WaveFunc.*WEIGHT_affix;
    Norb = length(orb_list);
    disp(Norb);
    %[Nwave,Nlist] = size(WaveFunc);
    if  Norb ~= Nwave
        error('Orbital list length is not equal to WaveFunc');
    end
    WFplot_list(Norb,:) = [0,0,0,0];
    for i = 1:Norb
        temp_C = 0;
        for j =1:Nlist
            temp_amp = WaveFunc(i,j);
            if temp_amp*temp_amp' < 1e-3
                temp_amp =0;
            end
            temp_C=temp_C+temp_amp*temp_amp';
        end
        WFamp = temp_C;
        WFplot_list(i,:) = [orb_list(i,:).*[1 1 1/10]*Rm ,WFamp];
    end
    figure();
    for j = 1:Ncycle
        temp_add = [0 0 (j-1)]*Rm;
        for i = 1 : Norb
            %plot(PS(j,1),PS(j,2),'ro','MarkerSize',(Z1(j)'*Z1(j) + Z2(j)'*Z2(j) + Z3(j)'*Z3(j) + Z4(j)'*Z4(j) + Z5(j)'*Z5(j) + Z6(j)'*Z6(j) +0.001)*200,'MarkerFaceColor','r');
            if WFplot_list(i,4) > 0
                
                plot3(WFplot_list(i,1)+temp_add(1),WFplot_list(i,2)+temp_add(2),WFplot_list(i,3)+temp_add(3),'ro','MarkerSize',(WFplot_list(i,4))*2000,'MarkerFaceColor','r');
            else
                plot3(WFplot_list(i,1)+temp_add(1),WFplot_list(i,2)+temp_add(2),WFplot_list(i,3)+temp_add(3),'ko','MarkerSize',1,'MarkerFaceColor','k');
            end
            hold on
        end
    end
    PARCHGCAR = WFplot_list;
    view(0,90);
    %%
elseif strcmp(mode,'PARCHG')
    cutoff = Ncycle;
    %% check
    disp(PARCHG_file);
    Norb = length(orb_list);
    disp(Norb);
    if  Norb ~= Nwave
        error('Orbital list length is not equal to WaveFunc');
    end
    %% init
    %     adjustparm = 30;
    %adjustparm = round(3/norm(curve_R))+1;
    
    a = norm(Rm(1,:));
    b = norm(Rm(2,:));
    c = norm(Rm(3,:));
    % Rmesh init
    NGXF = round(a/curve_R(1));
    NGYF = round(b/curve_R(2));
    NGZF = round(c/curve_R(3));
    for iz = 1:NGZF
        PARCHGCAR{iz} = sparse(NGXF,NGYF);
    end
    adjustparmX = round(cutoff/curve_R(1));
    adjustparmY = round(cutoff/curve_R(2));
    adjustparmZ = round(cutoff/curve_R(3));
    %     disp(NGXF);
    %     disp(NGYF);
    %     disp(NGZF);
    %% Data
    % init PARCHG_line
%     PARCHG_line_single = struct('seq',[],'NX',[],'NY',[],'NZ',[],'Rc',[],...
%         'To_orbital_list',[],'To_orbital_seq_list',[],'To_orbital_Rlength_list',[],'To_orbital_Wave_amp_list',[],'addornot',0);
%     PARCHG_line = repmat(PARCHG_line_single,[NGXF*NGYF*NGZF,1]);
    
    tic;
    for iorb = 1:Norb
        temp_orb = (orb_list(iorb,:));
        Ro(1) = round(temp_orb(1)*NGXF);
        Ro(2) = round(temp_orb(2)*NGYF);
        Ro(3) = round(temp_orb(3)*NGZF);
        %Ro = round(temp_orb./curve_R);
        %disp(Ro);
        
        for kz = -adjustparmZ : adjustparmZ
            for jy = -adjustparmY : adjustparmY
                for ix = -adjustparmX : adjustparmX
                    i = ix + Ro(1);
                    if i >0 && i<  NGXF
                        j = jy + Ro(2);
                        if j >0 &&j<  NGYF
                            k = kz + Ro(3);
                            if k>0 && k<  NGZF
                                Rc  = [(ix-1)/NGXF,(jy-1)/NGYF,(kz-1)/NGZF];
                                PARCHGCAR{k}(i,j) = PARCHGCAR{k}(i,j) + S_wave_amp(norm(Rc*Rm),WaveFunc(iorb,:)*10000,cutoff) ;
                            end
                            %disp(PARCHG_line(seq).Rc);
                        end
                    end
                end
            end
        end
        fprintf('(%6d/%-6d)th orb information has been added in PARCHGCAR.\n',iorb,Norb);
    end
    toc;
    fprintf('PARCHGCAR gen complete.\n');
    %% Write
    tic;
    fprintf(' Start write PARCHG.\n');
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
        temp_PARCHGCAR = full(PARCHGCAR{k});
        for j = 1:NGYF
            for i = 1:NGXF
                %seq = (k-1)*NGYF + (j-1)*NGXF + i;
                
                Ten_cont = Ten_cont + 1;
                %fprintf(fileID,"%8.4f  ", PARCHG_line(seq).C);
                %fprintf("%8.6f  ", temp_PARCHGCAR(i,j));
                fprintf(fileID,"%8.6f  ", temp_PARCHGCAR(i,j));
                if Ten_cont == 10
                    Ten_cont = 0;
                    fprintf(fileID,"\n");
                    
                end
            end
        end
        fprintf('The (%6d/%-6d) th PARCHG has been wirtten.\n',k,NGZF);
    end
    fclose(fileID);
    toc;
end
end


function C = S_wave_amp(Rr,Wave_amp,cutoff)
    % nargin
    % if nargin < 3
    %     adjustparm = 20;
    % end
    % check
    % if NWave_amp ~= Rlength_list
    %     error('nRlength_list ~= nWave_amp, why??');
    % end
        if Rr <  1
            Rr =1;
        end
        temp_C=0;
        for j =1:length(Wave_amp)
            temp_amp = Wave_amp(1,j);
%             if temp_amp*temp_amp' < 1e-4
%                 temp_amp =0exp
%             end
            temp_C=temp_C+temp_amp*temp_amp'*exp(-Rr/cutoff)^2;
        end
        C = temp_C;
end

