function [SG,R_struct] = cif_gen(filename,Rm,Atom_name,sites,SG,mode)
    if nargin < 1
        filename = "POSCAR.cif";
    end
    if nargin < 2
        [Rm,sites,Atom_name,~,~]=POSCAR_readin('POSCAR','vasp');
    end
    if nargin < 5
        SG = 1;
    end
    
    if nargin < 6
        mode = 'P1_cif';
    end    
    title = "# CIF Gen by MATLAB";
    fileID = fopen(filename,'w');
    fprintf(fileID,"%s\n",title);
    fprintf(fileID,'#email: parkman@buaa.edu.cn\n');
    fprintf(fileID,'## Date: %s\n\n',date);
    %
    fprintf(fileID,'# %s\n\n',filename);
    if strcmp(mode,'P1_mcif')
    % --- Displacement Data ---
        fprintf(fileID,'# --- Displacement Data ---\n\n');
        fprintf(fileID,"data_"+"park"+string(rand()*1000)+"\n");
        fprintf(fileID,'_audit_creation_date    %s\n'        ,   date);
        fprintf(fileID,'_audit_creation_method   %s\n'       ,   "vasplib");
        fprintf(fileID,'###non-st# _data_type "displacement"\n\n');
        % wanring
        fprintf(fileID,'#####################################################\n');
        fprintf(fileID,'# Please be warned that this structure is actually  #\n');
        fprintf(fileID,'# *NOT* magnetic but is represented as such to be   #\n');
        fprintf(fileID,'# able to visualize the symmetry modes.             #\n');
        fprintf(fileID,'#####################################################\n\n');
        fprintf(fileID,'_space_group.magn_number_BNS       "?"\n');
        fprintf(fileID,'_space_group.magn_name_BNS         "?"\n');
    end
    % CRYSTAL INFOR
    R_struct = Rm2abc(Rm);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_a",    R_struct.a);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_b",   R_struct.b);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_c", R_struct.c);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_angle_alpha", R_struct.alpha);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_angle_beta",  R_struct.beta);
    fprintf(fileID,'%-20s %13.6f\n\n',"_cell_angle_gamma ",R_struct.gamma);
    
    fprintf(fileID,'loop_\n');
    fprintf(fileID,'_space_group_symop.magn_id\n');
    fprintf(fileID,'_space_group_symop.magn_operation_xyz\n');
    fprintf(fileID,'1   x,y,z,+1\n\n');
    % SITES INFORMATION
    fprintf(fileID,'loop_\n');
    fprintf(fileID,'_atom_site_label\n');
    fprintf(fileID,'_atom_site_type_symbol\n');
    fprintf(fileID,'_atom_site_fract_x\n');
    fprintf(fileID,'_atom_site_fract_y\n');
    fprintf(fileID,'_atom_site_fract_z\n');
    fprintf(fileID,'_atom_site_occupancy\n');
    switch mode
        case {'P1_mcif','P1_cif'}
            for i =1:length(sites)
                fprintf(fileID,'%4s %3s %7.5f %7.5f %7.5f %6.4f\n',...
                    sites(i).name,Atom_name(sites(i).nameseq),...
                    sites(i).rc1,sites(i).rc2,sites(i).rc3,...
                    1.0000);
            end
        case {'CEF_cif'}
            for i =1:length(sites)
                fprintf(fileID,'%4s %3s %7.5f %7.5f %7.5f %6.4f\n',...
                    element2ion(Atom_name(sites(i).nameseq)),Atom_name(sites(i).nameseq),...
                    sites(i).rc1,sites(i).rc2,sites(i).rc3,...
                    1.0000);
            end            
            
    end
    fprintf(fileID,'\n');
    if strcmp(mode,'P1_mcif')
        fprintf(fileID,'loop_\n');
        fprintf(fileID,'_atom_site_moment_label\n');
        fprintf(fileID,'_atom_site_moment_crystalaxis_x\n');
        fprintf(fileID,'_atom_site_moment_crystalaxis_y\n');
        fprintf(fileID,'_atom_site_moment_crystalaxis_z\n');
        for i =1:length(sites)
            if sites(i).mag_ion
                if length(sites(i).mag_amp) == 3
                    fprintf(fileID,'%4s %7.5f %7.5f %7.5f\n',...
                        sites(i).name,...
                        R_struct.a*sites(i).mag_amp(1)/10,...
                        R_struct.b*sites(i).mag_amp(2)/10,...
                        R_struct.c*sites(i).mag_amp(3)/10);
                else
                    fprintf(fileID,'%4s %7.5f %7.5f %7.5f\n',...
                        sites(i).name,...
                        R_struct.a*sites(i).mag_amp(1)/10,...
                        R_struct.b*sites(i).mag_amp(1)/10,...
                        R_struct.c*sites(i).mag_amp(1)/10);
                end
            end

        end
    end
    fclose(fileID);
end
function R_struct = Rm2abc(Rm)
    R_struct.a = norm(Rm(1,:));
    R_struct.b = norm(Rm(2,:));
    R_struct.c = norm(Rm(3,:));
    R_struct.alpha = acos(dot(Rm(2,:),Rm(3,:))/(norm(R_struct.b)*norm(R_struct.c)))/pi*180;
    R_struct.beta  = acos(dot(Rm(3,:),Rm(1,:))/(norm(R_struct.c)*norm(R_struct.a)))/pi*180;
    R_struct.gamma = acos(dot(Rm(1,:),Rm(2,:))/(norm(R_struct.a)*norm(R_struct.b)))/pi*180;
end 
function ionname = element2ion(element)
    switch char(element)
        case {'Rb','Na','Li','K','H','Fr','Cs','Ag'} 
            ionname = strcat(element,'1+');
            
        case {'Be','Mg','Ca','Sr','Ba','Ra','Cu','Zn','Fe',...
            'Sm','Eu',...
            'Co','Ni','Ru','Pd','Pt','Hg','Cd'...
            } 
            ionname = strcat(element,'2+');
            
        case {'B','Al','Ga','In','Tl','Nh','La','Nd','Pm',...
            'Gd','Ho','Er','Tm','Yb','Lu',...   
            'Cr','Sc','Y','Rh','Ir','Au'...
            }
            ionname = strcat(element,'3+');
            
        case {'C','Si','Ge','Sn','Pb','Fl',...
            'Ce','Pr','Tb','Dy',...  
            'Ti','Zr','Hf','Rf','Mo','W','Sg','Os'...
            }
            ionname = strcat(element,'4+');
            
            
        case {'N','P',...    
                'V','Nb','Ta','Db'...
            }
            ionname = strcat(element,'5+');
            
        case {'Mn','Tc','Re','Bh'}
            ionname = strcat(element,'7+');
            
        case {'F','Cl','Br','I','At','Ts'}
            ionname = strcat(element,'1-') ;
            
        case {'O','S','Se','Te','Po','Lv'}
            ionname = strcat(element,'2-');    
            
        case {'As','Sb','Bi','Ms'}
            ionname = strcat(element,'3-');

        otherwise
              ionname = element;  
    end
        
end