load('Ising_TriAFM_exact_CvT.mat')

%% You can use another Cv T
% 给出数据
%   在工具栏给出 Cv T
Cv_origin = evalin('base','Cv');
T_origin  = evalin('base','T');
choose_data = 1:20;
Cv = Cv_origin(choose_data);
T = T_origin(choose_data);

%% 通过fit工具箱直接拟合
p1 = fittype('a*x^n');
p2 = fittype('c*exp(d*x)');
p3 = fittype('a*x^n+c*exp(d*x)');
disp(['代数拟合参数' ;coeffnames(p1)]);
disp(['指数拟合参数' ;coeffnames(p2)]);
disp(['指数-代数拟合参数'; coeffnames(p3)]);

f1 = fit(T,Cv,p1);
f2 = fit(T,Cv,p2);
f3 = fit(T,Cv,p3,'Lower',[0.0000,0,0,0],'Upper',[inf,inf,inf,5]);
disp('代数拟合结果');
disp(f1);
disp('指数拟合结果')
disp(f2);
disp('指数拟合结果')
disp(f3);
% 对比结果 永字八法
Cv1 = f1.a .* T.^ f1.n;
Cv2 = f2.c .* exp( f2.d*T);
Cv3 = f3.a .* T.^ f3.n + f3.c .* exp( f3.d*T);
% plot
fig = figure('PaperType','a4letter','PaperSize',[16 16],'Color','white');
% for
ax1 = subplot(2,2,1,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica");
ax2 = subplot(2,2,2,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica");
ax3 = subplot(2,2,3,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica");
ax4 = subplot(2,2,4,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica");
% end

box(ax1,'on');
hold(ax1,'all');
grid(ax1,'on')

box(ax2,'on');
hold(ax2,'all');
grid(ax2,'on')

box(ax3,'on');
hold(ax3,'all');
grid(ax3,'on')

box(ax4,'on');
hold(ax4,'all');
grid(ax4,'on');

semilogx(ax1,T,Cv,'DisplayName','C_v origin','LineWidth',2);
semilogx(ax1,T,Cv1,'DisplayName','代数: C_v = '+string(f1.a)+"*T^{"+string(f1.n)+"}",'LineWidth',2);
semilogx(ax1,T,Cv2,'DisplayName','指数: C_v = '+string(f2.c)+"*e^{"+string(f2.d)+'*T}','LineWidth',2);
semilogx(ax1,T,Cv3,'DisplayName','代数-指数混合: C_v = '+...
string(f3.a)+"*T^{"+string(f3.n)+"}+"+string(f3.c)+'*e^{'+string(f3.d)+'*T}','LineWidth',2);
legend1 = legend(ax1);
set(legend1,'FontSize',12);
set(ax1,'XScale','log') ;
ylabel(ax1,'C_v (J/K)');
title(ax1,'log -> T');

semilogy(ax2,T,Cv,'DisplayName','C_v origin','LineWidth',2);
semilogy(ax2,T,Cv1,'DisplayName','代数','LineWidth',2);
semilogy(ax2,T,Cv2,'DisplayName','指数','LineWidth',2);
%semilogy(ax2,T,Cv3,'DisplayName','代数-指数混合','LineWidth',2);
legend(ax2);
set(ax2,'YScale','log') ;
title(ax2,'log -> C_v');

loglog(ax3,T,Cv,'DisplayName','C_v origin','LineWidth',2);
loglog(ax3,T,Cv1,'DisplayName','代数','LineWidth',2);
loglog(ax3,T,Cv2,'DisplayName','指数','LineWidth',2);
%loglog(ax3,T,Cv3,'DisplayName','代数-指数混合','LineWidth',2);
legend(ax3);
set(ax3,'XScale','log') ;
set(ax3,'YScale','log') ;
ylabel(ax3,'C_v (J/K)');
xlabel(ax3,'T (K)');
title(ax3,'loglog')

plot(ax4,T.^-1,log(Cv),'DisplayName','C_v origin','LineWidth',2);
plot(ax4,T.^-1,log(Cv1),'DisplayName','代数','LineWidth',2);
plot(ax4,T.^-1,log(Cv2),'DisplayName','指数','LineWidth',2);
%plot(ax4,T.^-1,log(Cv3),'DisplayName','代数-指数混合','LineWidth',2);
legend(ax4);
xlabel(ax4,'T (K)');
title(ax4,'log(C_v) vs 1/T');