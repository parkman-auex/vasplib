function varargout =   OUTCAR_read(filename,mode)
    %
    arguments
        filename = "OUTCAR";
        mode = 'MAG';
    end
    import park.*;
    % SYSTEM
    [~,SYSTEM_str] = system("grep SYSTEM "+ filename);
    SYSTEM_str = strsplit(SYSTEM_str,'SYSTEM = ');
    SYSTEM_str = strsplit(strtrim(SYSTEM_str{2}));
    SYSTEM = string(SYSTEM_str{1});
    % NIONS
    [~,NIONS_str] = system("grep NIONS "+ filename);
    NIONS_str = strsplit(NIONS_str,'NIONS = ');
    NIONS = str2double(NIONS_str{2});
    % NBANDS
    [~,NBANDS_str] = system("grep NBANDS "+ filename);
    NBANDS_str = strsplit(NBANDS_str,'NBANDS');
    NBANDS_str = strsplit(NBANDS_str{2},'=');
    NBANDS = str2double(NBANDS_str{2});
    % ENERGY TOTEN
    [~,ENERGY_TOTEN_str] = system("grep 'free  energy   TOTEN' "+ filename);
    if ~strcmp(ENERGY_TOTEN_str,"")
        ENERGY_TOTEN_str = strsplit(ENERGY_TOTEN_str,'=');
        ENERGY_TOTEN_str = strsplit(strtrim(ENERGY_TOTEN_str{2}));
        ENERGY_TOTEN = str2double(ENERGY_TOTEN_str{1});
    else
        ENERGY_TOTEN = 0;
    end
    % NELECT
    [~,NELECT_str] = system("grep 'NELECT' "+ filename);
    NELECT_str = strsplit(NELECT_str,'=');
    NELECT_str = strsplit(strtrim(NELECT_str{2}));
    NELECT = str2double(NELECT_str{1});
    % LDAUU
    [~,LDAUU_str] = system("grep LDAUU "+ filename);
    LDAUU_str = strsplit(LDAUU_str,'LDAUU = ');
    LDAUU = str2double(strsplit(strtrim(LDAUU_str{2})));
    DAUU = LDAUU(LDAUU>0);
    if isempty(DAUU )
        DAUU  = 0;
    end
    % LNONCOLLINEAR
    [~,LNONCOLLINEAR_str] = system("grep LNONCOLLINEAR "+ filename);
    LNONCOLLINEAR_str = strsplit(LNONCOLLINEAR_str,'LNONCOLLINEAR = ');
    LNONCOLLINEAR_str = strsplit(strtrim(LNONCOLLINEAR_str{2}));
    LNONCOLLINEAR_str = LNONCOLLINEAR_str{1};
    if contains(LNONCOLLINEAR_str,'T')
        LNONCOLLINEAR = 1;
    else
        LNONCOLLINEAR = 0;
    end
    % LSORBIT
    [~,LSORBIT_str] = system("grep LSORBIT "+ filename);
    LSORBIT_str = strsplit(LSORBIT_str,'LSORBIT = ');
    LSORBIT_str = strsplit(strtrim(LSORBIT_str{2}));
    LSORBIT_str = LSORBIT_str{1};
    if contains(LSORBIT_str,'T')
        LSORBIT = 1;
    else
        LSORBIT = 0;
    end
    if strcmp(mode,'MAG')
        % disp
        fprintf('U(eV) NIONS NCOL NBANDS ENERGY_TOT(eV) SYSTEM\n');
        fprintf('%5d %5d %4d %6d %13.6f %s\n',DAUU(1),NIONS,LNONCOLLINEAR,NBANDS,ENERGY_TOTEN,SYSTEM);
        if LNONCOLLINEAR == 1
            [~,MAG_start_line_str] = system("grep -n magnetization "+ filename+...
                "|tail -n 3 |awk -F ':' '{print $1}'");
            MAG_start_line = str2double(strsplit(string(MAG_start_line_str),'\n'));
            MAG_start_line = MAG_start_line(1:3);
            MAG_end_line = MAG_start_line+NIONS+5;
            % x
            [~,MAG_str{1}] =  system("sed -n '"+string(MAG_start_line(1)+4) +...
                ","+string(MAG_end_line(1)-2)+"p' "+filename );
            % y
            [~,MAG_str{2}] =  system("sed -n '"+string(MAG_start_line(2)+4) +...
                ","+string(MAG_end_line(2)-2)+"p' "+filename );
            % z
            [~,MAG_str{3}] =  system("sed -n '"+string(MAG_start_line(3)+4) +...
                ","+string(MAG_end_line(3)-2)+"p' "+filename );
            
            %disp(MAG_str);
            if  ENERGY_TOTEN < 0
                [MAG(:,:,1),~] = str2double_mat(MAG_str{1});
                [MAG(:,:,2),~] = str2double_mat(MAG_str{2});
                [MAG(:,:,3),~] = str2double_mat(MAG_str{3});
            else
                MAG(:,:,1) = zeros(NIONS,5,1);
                MAG(:,:,2) = zeros(NIONS,5,1);
                MAG(:,:,3) = zeros(NIONS,5,1);
            end
            %disp(MAG);
            fprintf('ION MAG_x(uB) MAG_y(uB) MAG_z(uB)\n');
            MAG_final_list = zeros(NIONS,3);
            for i = 1:NIONS
                fprintf('%3d %9.4f %9.4f %9.4f\n',i,MAG(i,end,1),MAG(i,end,2),MAG(i,end,3));
                MAG_final_list(i,:) = [MAG(i,end,1),MAG(i,end,2),MAG(i,end,3)];
            end
            varargout{1} = MAG_final_list;
            varargout{2} = DAUU;
            varargout{3} = ENERGY_TOTEN;
            varargout{4} = SYSTEM;
            varargout{5} = NELECT;
            varargout{6} = LSORBIT;
        else
            
            [~,MAG_start_line_str] = system("grep -n magnetization "+ filename+...
                "|tail -n 1 |awk -F ':' '{print $1}'");
            MAG_start_line = str2double(strsplit(string(MAG_start_line_str),'\n'));
            MAG_start_line = MAG_start_line(1);
            MAG_end_line = MAG_start_line+NIONS+5;
            % x
            [~,MAG_str] =  system("sed -n '"+string(MAG_start_line(1)+4) +...
                ","+string(MAG_end_line(1)-2)+"p' "+filename );
            if  ENERGY_TOTEN < 0
                [MAG(:,:,1),~] = str2double_mat(MAG_str);
            else
                MAG(:,:,1) = zeros(NIONS,5,1);
            end
            
            fprintf('ION MAG(uB) \n');
            MAG_final_list = zeros(NIONS,1);
            for i = 1:NIONS
                fprintf('%3d %9.4f\n',i,MAG(i,end,1));
                 MAG_final_list(i,:) = MAG(i,end,1);
            end
    
            varargout{1} = MAG_final_list;
            varargout{2} = DAUU;
            varargout{3} = ENERGY_TOTEN;
            varargout{4} = SYSTEM;
            varargout{5} = NELECT;
            varargout{6} = LSORBIT;
        end
    end
end

