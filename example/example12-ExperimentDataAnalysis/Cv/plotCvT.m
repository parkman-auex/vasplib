import vasplib_experiment.*;
[fig,ax] = creat_figure();
plot(ax,T,Cv,'LineWidth',3,'Color','b','DisplayName','origin');
plot(ax,T,Cv_prime,'LineWidth',3,'Color','r','DisplayName','fitting');
Cv_rm_N = Cv - Nulear_coe./(T.^2);
plot(ax,T,Cv_rm_N,'LineWidth',3,'Color','g','DisplayName','rm-C_v^{nuclear}');
legend();