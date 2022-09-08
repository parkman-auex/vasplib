function [WindingNumber,fig] = TopoCharge(H_hk,GammaOper,kloop,kpoints_r)
            if nargin < 3
                kloop = H_hk.kloop_gen();
            end
            if nargin < 2
                GammaOper = diag([ones(1,H_hk.Basis_num/2),-ones(1,H_hk.Basis_num/2)]);
            end
            if nargin < 4
                mode = 'norm';
            else
                mode = 'plot';
            end
            % 
            syms('k_x','k_y','k_z')
            TargetH =  H_hk.Hk_num ;
            TargetH_Pk_x = diff(TargetH,k_x);
            TargetH_Pk_y = diff(TargetH,k_y);
            TargetH_Pk_z = diff(TargetH,k_z);
            TargetHfun = matlabFunction(TargetH,'Vars',[k_x k_y k_z]);
            TargetHfun_Pk_x = matlabFunction(TargetH_Pk_x,'Vars',[k_x k_y k_z]);
            TargetHfun_Pk_y = matlabFunction(TargetH_Pk_y,'Vars',[k_x k_y k_z]);
            TargetHfun_Pk_z = matlabFunction(TargetH_Pk_z,'Vars',[k_x k_y k_z]);
            % Topo = intral(Gamme H-1 PH dk)
            kn = size(kloop,1);
            % uncheck loop
            dk_list = kloop(2:kn,:)- kloop(1:kn-1,:);
            dk_list = [dk_list;dk_list(1,:)];
            WindingNumber = 0;
            fig = figure('PaperType','a4letter','PaperSize',[16 8],'Color','white');
            % for

            ax0 = subplot(1,3,1,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica");
            box(ax0,'on');
            hold(ax0,'all');
            xlabel(ax0,'k_1');
            ylabel(ax0,'k_2');
            zlabel(ax0,'k_3');
            if strcmp(mode,'plot')
                h = kpoints_r(3); % 高度
                r = sqrt(kpoints_r(1)^2+kpoints_r(2)^2);  %半径
                if r == 0
                    plot3([0 0],[0 0],[-h h],'DisplayName','Node line','LineWidth',3);
                else
                    pos = [0,0,h]; % 圆心位置
                    t=0:0.001:(2*pi);  % 圆滑性设置
                    t=[t,0];
                    plot3(ax0,pos(1)+r*sin(t),pos(2)+r*cos(t), h*ones(size(t)),'DisplayName','Node line','LineWidth',3);
                end
               
            end
            title(ax0,'integral path');
            ax1 = subplot(1,3,2,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica");
            hold(ax1,'all');
            xlabel(ax1,'kpath num');
            ylabel(ax1,'dTopoCharge');
            ax2 = subplot(1,3,3,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica");
            hold(ax2,'all');
            xlabel(ax2,'kpath num');
            ylabel(ax2,'WindingNumber');
            dTopoCharge_old = 0;
            for i =1:kn
                kx = kloop(i,1);
                ky = kloop(i,2);
                kz = kloop(i,3);
                Ak = ...
                    +TargetHfun_Pk_x(kx,ky,kz)*dk_list(i,1)+...
                    +TargetHfun_Pk_y(kx,ky,kz)*dk_list(i,2)+...
                    +TargetHfun_Pk_z(kx,ky,kz)*dk_list(i,3);...
                    
                    Htemp = GammaOper/(TargetHfun(kx,ky,kz))*Ak;
                disp(Htemp);
                dTopoCharge = imag(trace(Htemp));
                WindingNumber =WindingNumber + dTopoCharge/(4*pi);
                scatter3(ax0,kx,ky,kz);
                area(ax1,[i-1,i],[dTopoCharge_old,dTopoCharge]);
                dTopoCharge_old = dTopoCharge;
                %scatter(ax1,i,dTopoCharge);
                stem(ax2,i,WindingNumber);
                %drawnow;
            end
            title(ax2,'WindingNumber = '+string(round(WindingNumber)));
            axis(ax0,'equal');
            %WindingNumber = WindingNumber;
        end