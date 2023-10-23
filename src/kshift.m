function [klist_new,DATA_L_new] = kshift(klist,KCUBE_BULK,DATA_L,opt)
arguments
    klist =[]; 
    KCUBE_BULK  = [-0.5 -0.5 -0.5;1 0 0;0 1 0;0 0 1];
    DATA_L = [];
    opt.cart = false;
    opt.Rm = POSCAR_read;
    opt.tol = 1e-12;
    opt.MeshOutput = false;
    opt.method {mustBeMember(opt.method,{'linear','nearest','natural'})} = 'linear';
end
%% DATA format
sizemesh = size(DATA_L);
if length(sizemesh) == 3
    DATA_L = reshape(DATA_L,[sizemesh(1)*sizemesh(2),sizemesh(3)]);
elseif ( (size(klist,1)*size(klist,2) == sizemesh(1)*sizemesh(2))|| (size(klist,1) == sizemesh(1)*sizemesh(2)) ) && length(sizemesh) == 2
    DATA_L = reshape(DATA_L,[sizemesh(1)*sizemesh(2),1]);
    %DATA_L = DATA_L(:);
end
%figure()
%scatter3(klist(:,1),klist(:,2),klist(:,3),abs(DATA_L)*10000,DATA_L);view(2);axis equal
%surf(Grid(:,:,1),Grid(:,:,2),Grid(:,:,3),reshape(DATA_L,[100 100]))
sizemeshk = size(klist);
if length(sizemeshk) == 3
    Grid = klist;
    klist = zeros(sizemesh(1)*sizemesh(2),3);
    klist(:,1) = reshape(Grid(:,:,1),[sizemesh(1)*sizemesh(2),1]);
    klist(:,2) = reshape(Grid(:,:,2),[sizemesh(1)*sizemesh(2),1]);
    klist(:,3) = reshape(Grid(:,:,3),[sizemesh(1)*sizemesh(2),1]);
    %Gridtmp = Grid(:,:,1);klist(:,1) = Gridtmp(:).';
    %Gridtmp = Grid(:,:,2);klist(:,2) = Gridtmp(:).';
    %Gridtmp = Grid(:,:,3);klist(:,3) = Gridtmp(:).';
    %figure()
    %surf(Grid(:,:,1),Grid(:,:,2),Grid(:,:,3),reshape(DATA_L,[100 100]))
end
%%
if opt.cart
    Gk = (2*pi*eye(3)/opt.Rm).';
    klist = klist/Gk;
else
    Gk = (2*pi*eye(3)/opt.Rm).';
end
%% unique!
[klist,uniqueseq,~] = uniquetol(klist,opt.tol,'ByRows',true);
if ~isempty(DATA_L)
    DATA_L = DATA_L(uniqueseq,:);
end
% figure()
% scatter3(klist(:,1),klist(:,2),klist(:,3),abs(DATA_L)*10000,DATA_L);view(2);axis equal
%% Converge KCUBE_BULK
switch size(KCUBE_BULK,1)
    case {0,1}
    case 2
    case 3 % 2D we use cart
        KCUBE_BULK(4,:) = KCUBE_BULK(2,:) + KCUBE_BULK(3,:);
        % Project mat too complex 
        % suppse pure 2D
        %K1 = KCUBE_BULK(2,:);K2 = KCUBE_BULK(3,:);
        %K1 = KCUBE_BULK(2,:);K2 = KCUBE_BULK(3,:);
        %Pplane_cart = [K1/norm(K1);K2/norm(K2)].';
        KCUBE_BULK(2:4,:) = KCUBE_BULK(2:4,:) + KCUBE_BULK(1,:);
        CENTER = (KCUBE_BULK(1,:) + KCUBE_BULK(4,:))/2;
        CENTER_base = round(CENTER); %remove  to origin;
        %CENTER_base_cart = CENTER_base*Gk;
        %KCUBE_BULK_cart = KCUBE_BULK*Gk;
        %KCUBE_BULK_project = KCUBE_BULK_cart*Pplane_cart;
        KCUBE_BULK = KCUBE_BULK(:,1:2);
        %klist_project = klist*Pplane;
        mode = '2D';
    case 4
        KCUBE_BULK(5,:) = KCUBE_BULK(2,:) + KCUBE_BULK(3,:);
        KCUBE_BULK(6,:) = KCUBE_BULK(2,:) + KCUBE_BULK(4,:);
        KCUBE_BULK(7,:) = KCUBE_BULK(3,:) + KCUBE_BULK(4,:);
        KCUBE_BULK(8,:) = KCUBE_BULK(2,:) + KCUBE_BULK(3,:) + KCUBE_BULK(4,:);
        KCUBE_BULK(2:8,:) = KCUBE_BULK(2:8,:) + KCUBE_BULK(1,:);
        CENTER = (KCUBE_BULK(1,:) + KCUBE_BULK(8,:))/2;
        CENTER_base = round(CENTER); %remove  to origin
        mode = '3D';
    otherwise
        %KCUBE_BULK(2:end,:) = KCUBE_BULK(2:end,:)+KCUBE_BULK(1,:);
end
%% creat a polyhedron
% https://www.mathworks.com/help/matlab/ref/delaunaytriangulation.convexhull.html?searchHighlight=convexhull&s_tid=srchtitle_convexhull_1
DT = delaunayTriangulation(KCUBE_BULK);% suppse pure 2D
[~,v] = convexHull(DT);
% Map ALL Klabel TO [0,1)  ?seems this procedure is totally wrong
klist = mod(klist,1);
% klist(klist == 1) = 0;   ?seems we dont need 1 -> 0
%Display the volume and plot the convex hull.
% trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3), ...
%        'FaceColor','cyan')
% auto V
Vceil = ceil(abs(v)); VceilL = (-Vceil:Vceil).';nVceilL = length(VceilL);
switch mode
    case '2D'
        SuperL= zeros(nVceilL^2,3);
        count = 0;
        for i = -Vceil:Vceil
            for j = -Vceil:Vceil
                count = count + 1;
                SuperL(count,:) = [i,j,0]-CENTER_base;
            end
        end
        % expand klist
        nklist = size(klist,1);
        SuperL = kron(SuperL,ones([nklist 1]));
        klist_new = repmat(klist,[count,1])+SuperL;
        %nChooseL = repmat((1:nklist).',[count,1]);
        %% How to check if a given point lies inside a polygon
        % https://www.mathworks.com/help/matlab/ref/triangulation.pointlocation.html
        ID = pointLocation(DT,klist_new(:,1:2));
        %klist_new = klist_new / Pplane;
        %klist = klist/Pplane;
    case '3D'
        SuperL= zeros(nVceilL^3,3);
        count = 0;
        for i = -Vceil:Vceil
            for j = -Vceil:Vceil
                for k = -Vceil:Vceil
                    count = count + 1;
                    SuperL(count,:) = [i,j,k]-CENTER_base;
                end
            end
        end
        %klist = klist - [0.5,0.5,0.5];
        % expand klist
        nklist = size(klist,1);
        SuperL = kron(SuperL,ones([nklist 1]));
        klist_new = repmat(klist,[count,1])+SuperL;
        %nChooseL = repmat((1:nklist).',[count,1]);
        %% How to check if a given point lies inside a polyhedron
        % https://www.mathworks.com/help/matlab/ref/triangulation.pointlocation.html
        ID = pointLocation(DT,klist_new);
end
% Map 
%nChooseL = nChooseL(~isnan(ID));
klist_new = klist_new(~isnan(ID),:);% +CENTER_base never shift back

klist_new_map = mod(klist_new,1);
% klist_new_map(klist_new_map == 1) = 0; ?seems we dont need 1 -> 0


[~,reseqL] = ismembertol(klist_new_map,klist,opt.tol,'ByRows',true);
rmlist = reseqL == 0;

if ~sum(rmlist)
    if ~isempty(DATA_L)
        DATA_L_new = DATA_L(reseqL,:);
%         figure()
%         scatter3(klist_new(:,1),klist_new(:,2),klist_new(:,3),abs(DATA_L_new)*10000,DATA_L_new);view(2);axis equal
    else
        DATA_L_new = DATA_L;
    end
else
    if isempty(DATA_L)
        DATA_L = ones(size(klist,1),1);
    end
    if ~opt.MeshOutput
        switch mode
            case '2D'
                % interp2!!
                DATA_L_new = zeros(size(klist_new_map,1),size(DATA_L,2));
                klist1 = klist(:,1); klist2 = klist(:,2);%klist3 = klist(:,3);
                klistnew1 = klist_new_map(:,1); klistnew2 = klist_new_map(:,2);%klistnew3 = klist_new_map(:,3);
                for i = 1:size(DATA_L,2)
                    F = scatteredInterpolant(klist1,klist2,DATA_L(:,i),opt.method);
                    %meshU = griddata(meshX,meshY,'v4');
                    DATA_L_new(:,i) = F(klistnew1,klistnew2);
                end
            case '3D'
                % interp3!!
                DATA_L_new = zeros(size(klist_new_map,1),size(DATA_L,2));
                klist1 = klist(:,1); klist2 = klist(:,2);klist3 = klist(:,3);
                klistnew1 = klist_new_map(:,1); klistnew2 = klist_new_map(:,2);klistnew3 = klist_new_map(:,3);
                for i = 1:size(DATA_L,2)
                    F = scatteredInterpolant(klist1,klist2,klist3,DATA_L(:,i));
                    %meshU = griddata(meshX,meshY,meshZ,'linear');
                    DATA_L_new(:,i) = F(klistnew1,klistnew2,klistnew3);
                end
        end
    else

    end
end
%
if opt.cart
    klist_new = klist_new*Gk;
end
end