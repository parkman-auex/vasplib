%% To generate k-mesh for DOSplot or 2D-Bandplot
% usage [klist_dos]=kmesh3D(Rm,mesh)
% input  : mesh [meshx meshy meshz]
%          Rm
            %a1 a1x a1y a1z
            %a2 a2x a2y a2z
            %a3 a3x a3y a3z 
% output : mesh [meshx meshy meshz]




function [klist_dos,klist1,klist2]=kmesh3D(Rm,mesh,kz,mode)
%--------  nargin  --------
if nargin <3
    mode = 'bulk-DOS';
    disp('kmesh 3D');
    disp('function needs to be repairaed')
elseif nargin <4
    mode = '3Dbandplot';
    disp('3D band plot');
end
%--------  init  --------
Gk  = (eye(3)*2*pi/Rm)' ;% Reciprocal vector
%--------  distri  --------
switch mode
    case 'bulk-DOS'
        b_x=sum(Gk(:,1));
        b_y=sum(Gk(:,2));
        b_z=sum(Gk(:,3));
        xnodes=mesh(1)*2+1;
        ynodes=mesh(2)*2+1;
        znodes=mesh(3)*2+1;
        xlist=linspace(-b_x,b_x,xnodes);
        ylist=linspace(-b_y,b_y,ynodes);
        zlist=linspace(-b_z,b_z,znodes);
        klist_dos=[];
        for i=1:xnodes
            for j=1:ynodes
                for k=1:znodes
                    klist_temp=[xlist(i) ylist(j) zlist(k)];
                    klist_dos=[klist_dos;klist_temp];
                end
            end
        end
        klist1 = xlist;
        klist2 = ylist;
    case '3Dbandplot'
        disp('Primitive cell')
        disp(Rm);
        disp('Reciprocal Lattice')
        disp(Gk);
        disp("Rm*Gk' ");
        disp(sym(Rm*Gk'));
        xnodes=mesh(1)*2+1;
        ynodes=mesh(2)*2+1;                
        klist1=linspace(-0.5,0.5,xnodes);
        klist2=linspace(-0.5,0.5,ynodes);

        %     kz_r = kz*Gk(3,:);
        klist_dos=[];
        for i=1:xnodes
            for j=1:ynodes
                klist_temp=[klist1(i) klist2(j) kz]*Gk;
                klist_dos=[klist_dos;klist_temp];
            end
        end
    case 'fermiarc'
        disp('fermi-arc');
        xnodes = mesh(1);
        ynodes = mesh(2);
        kfermiarc = kz;
        klist1 = [linspace(0,kfermiarc(2,1),xnodes)',...
                  linspace(0,kfermiarc(2,2),xnodes)',...
                  linspace(0,kfermiarc(2,3),xnodes)',...
                  ];
        klist2 = [linspace(0,kfermiarc(3,1),ynodes)',...
                  linspace(0,kfermiarc(3,2),ynodes)',...
                  linspace(0,kfermiarc(3,3),ynodes)',...
                  ];
        klist_1 = klist1 ;
        klist_2 = klist2 ;
%         disp( klist_1 );
%         disp( klist_2 );
        klist_dos = zeros(xnodes*ynodes,3);
        count =0;
        for i=1:xnodes
            for j=1:ynodes
                
                count = count+1;
                klist_dos(count,:) = klist_1(i,:)+klist_2(j,:)+ kfermiarc(1,:);
            end
        end
%         klist1 = linspace(kfermiarc(1,1),max(kfermiarc(2,:)),xnodes)';
%         klist2 = linspace(kfermiarc(1,2),max(kfermiarc(3,:)),xnodes)';
end

end






