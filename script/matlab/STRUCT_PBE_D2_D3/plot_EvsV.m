[Fig,ax] = vasplib_tool.create_figure(1);
ax1 = ax(1);
%%
load('EvsV_PBE.mat')
%  Fig = Figs(1,2)
%  ax1 = Fig.axes(1)
%  ax2 = Fig.axes(2)
set(ax1,'FontName','Times New Roman');
% set(ax2,'FontName','Times New Roman');
ylabel(ax1,'Energy (eV/atom)','interpreter','latex');
% ylabel(ax2,'Pressure (Gpa)','interpreter','latex');
xlabel(ax1,'Volume ({\AA}$^3$/atom)','interpreter','latex');
% xlabel(ax2,'Volume ({\AA}$^3$/atom)','interpreter','latex');
Func = 'PBE';
scatter(ax1,TOTVlist_r,TOTElist,100,'o','filled','DisplayName',"DFT-"+Func,'Marker','o');
%%
load('EvsV_D2.mat')
%  Fig = Figs(1,2)
%  ax1 = Fig.axes(1)
%  ax2 = Fig.axes(2)

set(ax1,'FontName','Times New Roman');
% set(ax2,'FontName','Times New Roman');
ylabel(ax1,'Energy (eV/atom)','interpreter','latex');
% ylabel(ax2,'Pressure (Gpa)','interpreter','latex');
xlabel(ax1,'Volume ({\AA}$^3$/atom)','interpreter','latex');
% xlabel(ax2,'Volume ({\AA}$^3$/atom)','interpreter','latex');
Func = 'D2';
scatter(ax1,TOTVlist_r,TOTElist,100,'o','filled','DisplayName',"DFT-"+Func,'Marker','v');
%%
load('EvsV_D3.mat')
%  Fig = Figs(1,2)
%  ax1 = Fig.axes(1)
%  ax2 = Fig.axes(2)

set(ax1,'FontName','Times New Roman');
% set(ax2,'FontName','Times New Roman');
ylabel(ax1,'Energy (eV/atom)','interpreter','latex');
% ylabel(ax2,'Pressure (Gpa)','interpreter','latex');
xlabel(ax1,'Volume ({\AA}$^3$/atom)','interpreter','latex');
% xlabel(ax2,'Volume ({\AA}$^3$/atom)','interpreter','latex');
Func = 'D3';
scatter(ax1,TOTVlist_r,TOTElist,100,'o','filled','DisplayName',"DFT-"+Func,'Marker','square');
%%
legend();
export_fig EvsV.png -r300
save('plotresult.mat','Fig','ax','ax1');