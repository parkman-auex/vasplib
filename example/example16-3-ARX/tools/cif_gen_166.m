function R_struct = cif_gen_166(filename,Rm,Atom_name,sites,mode)
    if nargin < 1
        filename = "POSCAR.cif";
    end
    if nargin < 2
        [Rm,sites,Atom_name,~,~]=POSCAR_readin('POSCAR','vasp');
    end
    if nargin < 5
        mode = 'P1_cif';
    end
    title = "# CIF Gen by MATLAB";
    fileID = fopen(filename,'w');
    fprintf(fileID,"%s\n",title);
    fprintf(fileID,'#email: parkman@buaa.edu.cn\n');
    fprintf(fileID,'## Date: %s\n\n',date);
    %
    fprintf(fileID,'# %s\n\n',filename);
    
    R_struct = Rm2abc_166(Rm);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_a",    R_struct.a);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_b",   R_struct.b);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_c", R_struct.c);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_angle_alpha", R_struct.alpha);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_angle_beta",  R_struct.beta);
    fprintf(fileID,'%-20s %13.6f\n\n',"_cell_angle_gamma ",R_struct.gamma);
    fprintf(fileID,'%-20s %13.6f\n\n',"_cell_volume ",R_struct.volume);
    fprintf(fileID,'_space_group_name_H-M_alt              ''R -3 m'' \n');
    fprintf(fileID,'_space_group_IT_number                 166 \n');  
fprintf(fileID,'\n');    
fprintf(fileID,'loop_ \n');  
fprintf(fileID,'_space_group_symop_operation_xyz \n');  
fprintf(fileID,"   'x, y, z' \n");
fprintf(fileID,"   '-x, -y, -z'\n");
fprintf(fileID,"   '-y, x-y, z'\n");
fprintf(fileID,"   'y, -x+y, -z'\n");
fprintf(fileID,"   '-x+y, -x, z'\n");
fprintf(fileID,"   'x-y, x, -z'\n");
fprintf(fileID,"   'y, x, -z'\n");
fprintf(fileID,"   '-y, -x, z'\n");
fprintf(fileID,"   'x-y, -y, -z'\n");
fprintf(fileID,"   '-x+y, y, z'\n");
fprintf(fileID,"   '-x, -x+y, -z'\n");
fprintf(fileID,"   'x, x-y, z'\n");
fprintf(fileID,"   'x+2/3, y+1/3, z+1/3'\n");
fprintf(fileID,"   '-x+2/3, -y+1/3, -z+1/3'\n");
fprintf(fileID,"   '-y+2/3, x-y+1/3, z+1/3'\n");
fprintf(fileID,"   'y+2/3, -x+y+1/3, -z+1/3'\n");
fprintf(fileID,"   '-x+y+2/3, -x+1/3, z+1/3'\n");
fprintf(fileID,"   'x-y+2/3, x+1/3, -z+1/3'\n");
fprintf(fileID,"   'y+2/3, x+1/3, -z+1/3'\n");
fprintf(fileID,"   '-y+2/3, -x+1/3, z+1/3'\n");
fprintf(fileID,"   'x-y+2/3, -y+1/3, -z+1/3'\n");
fprintf(fileID,"   '-x+y+2/3, y+1/3, z+1/3'\n");
fprintf(fileID,"   '-x+2/3, -x+y+1/3, -z+1/3'\n");
fprintf(fileID,"   'x+2/3, x-y+1/3, z+1/3'\n");
fprintf(fileID,"   'x+1/3, y+2/3, z+2/3'\n");
fprintf(fileID,"   '-x+1/3, -y+2/3, -z+2/3'\n");
fprintf(fileID,"   '-y+1/3, x-y+2/3, z+2/3'\n");
fprintf(fileID,"  'y+1/3, -x+y+2/3, -z+2/3'\n");
fprintf(fileID,"   '-x+y+1/3, -x+2/3, z+2/3'\n");
fprintf(fileID,"   'x-y+1/3, x+2/3, -z+2/3'\n");
fprintf(fileID,"   'y+1/3, x+2/3, -z+2/3'\n");
fprintf(fileID,"   '-y+1/3, -x+2/3, z+2/3'\n");
fprintf(fileID,"   'x-y+1/3, -y+2/3, -z+2/3'\n");
fprintf(fileID,"   '-x+y+1/3, y+2/3, z+2/3'\n");
fprintf(fileID,"   '-x+1/3, -x+y+2/3, -z+2/3'\n");
fprintf(fileID,"   'x+1/3, x-y+2/3, z+2/3'\n");
fprintf(fileID,'\n');  
    % SITES INFORMATION 
       
    fprintf(fileID,'loop_\n');
    fprintf(fileID,'_atom_site_label\n');
    fprintf(fileID,'_atom_site_occupancy\n');
    fprintf(fileID,'_atom_site_fract_x\n');
    fprintf(fileID,'_atom_site_fract_y\n');
    fprintf(fileID,'_atom_site_fract_z\n');
     fprintf(fileID,' _atom_site_adp_type\n');
    fprintf(fileID,'_atom_site_U_iso_or_equiv\n');
    fprintf(fileID,'_atom_site_type_symbol\n');
    switch mode
        case {'P1_mcif','P1_cif'}
            for i =1:length(sites)
                fprintf(fileID,'%4s %3s %7.5f %7.5f %7.5f %6.4f\n',...
                    sites(i).name,Atom_name(sites(i).nameseq),...
                    sites(i).rc1,sites(i).rc2,sites(i).rc3,...
                    1.0000);
            end
        case {'CEF_cif'}
            for i =1:length(sites)-1
                fprintf(fileID,'%3s %3.1f %7.5f %7.5f %7.5f %s %s %4s\n',...
                    Atom_name(sites(i).nameseq),1,...
                    0,0,sites(i).rc3,...
                    'Uiso' , '?', ...
                    element2ion(Atom_name(sites(i).nameseq))...
                    );
            end            
            
    end
    
end
function R_struct = Rm2abc_166(Rm)
    Ns = [1 -1 0;0 1 -1;1 1 1];
    Rm = Ns*Rm;
    R_struct.a = norm(Rm(1,:));
    R_struct.b = norm(Rm(2,:));
    R_struct.c = norm(Rm(3,:));
    R_struct.alpha = acos(dot(Rm(2,:),Rm(3,:))/(norm(R_struct.b)*norm(R_struct.c)))/pi*180;
    R_struct.beta  = acos(dot(Rm(3,:),Rm(1,:))/(norm(R_struct.c)*norm(R_struct.a)))/pi*180;
    R_struct.gamma = acos(dot(Rm(1,:),Rm(2,:))/(norm(R_struct.a)*norm(R_struct.b)))/pi*180;
    R_struct.volume = norm(Rm);
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