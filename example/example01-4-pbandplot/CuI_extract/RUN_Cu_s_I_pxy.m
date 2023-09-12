%RUN_Cu_s_I_pxy
import vasplib_tool.*;
[WEIGHTCAR_struct_cell,Name_list,~] = WEIGHTCAR_gen_dir('PBAND');
WEIGHTCAR_struct_Cu = WEIGHTCAR_struct_cell{1};
WEIGHTCAR_struct_I1 = WEIGHTCAR_struct_cell{2};
WEIGHTCAR_struct_I2 = WEIGHTCAR_struct_cell{3};
WEIGHTCAR_struct_plot(1) =  WEIGHTCAR_struct_Cu(1);
WEIGHTCAR_struct_plot(2) =  WEIGHTCAR_struct_I1(2); %py
WEIGHTCAR_struct_plot(3) =  WEIGHTCAR_struct_I1(4); %px
% plot
EIGENCAR = EIGENVAL_read;
Ecut = [-3,3];
titlestring = 'fatband-CuI-nosoc';
% [fig,ax]=pbandplot(WEIGHTCAR_struct_plot,EIGENCAR,Ecut,titlestring);

% for another plot
WEIGHTCAR_struct_plot_2(1) =  WEIGHTCAR_struct_Cu(1);
WEIGHTCAR_struct_plot_2(2) =  WEIGHTCAR_struct_I1(2); 
WEIGHTCAR_struct_plot_2(2).displayname =  'I1-px,py'; %pxy
WEIGHTCAR_struct_plot_2(2).WEIGHTCAR = WEIGHTCAR_struct_I1(2).WEIGHTCAR+ ...
                                       WEIGHTCAR_struct_I1(4).WEIGHTCAR;%pxy
WEIGHTCAR_struct_plot_2(2).WEIGHTCAR =  WEIGHTCAR_struct_plot_2(2).WEIGHTCAR/2; 
% adjust density
%% 
[fig,ax]=vasplib_plot.pbandplot(WEIGHTCAR_struct_plot_2,EIGENCAR,'Ecut',Ecut,'title',titlestring);     

AxChildren = ax.Children
legend(AxChildren([4,5]))