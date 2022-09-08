%% An example for the visualization of 3D-bandstructure for your own by using MATLAB
% usage: figure_k = bandplot_3d(bandlist,Ecut,contourf_mode,wt_mode)
function figure_k = bandplot_3d(bandlist,Ecut,kplane,contourf_mode,wt_mode)
figure_k = figure();
if nargin <1
    bandlist = [1,2];
end
if nargin<2
    Ecut = [-1.5,1.5];
end
if nargin <3
    kplane = [1,2];
end
if nargin<4
    contourf_mode = 0;
end
if nargin<5
    wt_mode = 0;
end
% Create axes
axes1 = axes('Parent',figure_k );
hold(axes1,'on');
set(axes1,'FontName','Helvetica','FontSize',30);

if wt_mode == 0
    kx_mesh=load('KX.grd');      
    ky_mesh=load('KY.grd');   

    num_bandlist = length(bandlist);
    count_num = 0;
    BM_mesh{count_num+1} = ky_mesh;
    for i = bandlist
    count_num = count_num+1;
    temp_name = strcat('BAND.B',num2str(i),'.grd');
    BM_mesh{count_num}=load(temp_name);  
    end 

    for i = 1:num_bandlist

    surf(kx_mesh,ky_mesh,BM_mesh{i},'Parent',axes1,'FaceColor','interp',...
        'EdgeColor','none','FaceAlpha','0.6');
    hold on

    end

    if contourf_mode ==1
    for i = 1:num_bandlist

    contourf(kx_mesh,ky_mesh,BM_mesh{i},4)
    hold on

    end
    end
    xlabel('$\it{k}_{x} (\textrm{1/\AA})$','Interpreter','latex','Fontname','Times New Roman')
    ylabel('$\it{k}_{y} (\textrm{1/\AA})$','Interpreter','latex','Fontname','Times New Roman')
    zlabel('Energy (eV)')
    axis image
    axis vis3d
    shading interp;
    colormap(hsv);
    set(gca,'FontSize',20);
    zlim(Ecut);
elseif wt_mode == 1
    %% wt mode
    bulkek_plane_matlab = load('bulkek_plane-matlab.dat');
    [~,width]= size(bulkek_plane_matlab );
    kx_list = bulkek_plane_matlab(:,1);
    ky_list = bulkek_plane_matlab(:,2);
    kz_list = bulkek_plane_matlab(:,3);
    k1_list = bulkek_plane_matlab(:,4);
    k2_list = bulkek_plane_matlab(:,5);
    k3_list = bulkek_plane_matlab(:,6);
    num_k = sqrt(length(k2_list));
    
    BM_mesh{1} = reshape(bulkek_plane_matlab(:,7),num_k,num_k);
    BM_mesh{2} = reshape(bulkek_plane_matlab(:,8),num_k,num_k);
    if width >8
    BM_mesh{3} = reshape(bulkek_plane_matlab(:,9),num_k,num_k);
    BM_mesh{4} = reshape(bulkek_plane_matlab(:,10),num_k,num_k);
    end
    
    kx_mesh = reshape(bulkek_plane_matlab(:,kplane(1)),num_k,num_k);
    ky_mesh = reshape(bulkek_plane_matlab(:,kplane(2)),num_k,num_k);
    
    
    for i = bandlist

    surf(kx_mesh,ky_mesh,BM_mesh{i},'Parent',axes1,'FaceColor','interp',...
        'EdgeColor','none','FaceAlpha','0.6');
    hold on

    end
    
    if contourf_mode ==1
    for i = bandlist

    contourf(kx_mesh,ky_mesh,BM_mesh{i},4)
    hold on

    end
    end

    xlabel('k_x ','Fontname','Helvetica')
    ylabel('k_y','Fontname','Helvetica')
    zlabel('Energy (eV)')
    axis image
    axis vis3d
    shading interp;
    colormap(hsv);
    set(gca,'FontSize',30);
    zlim(Ecut);
    view(axes1,[-72 17]);
    % Create title
    dirname=pwd;
    dirname=strsplit(dirname,'/');
    titlestring=dirname(length(dirname));
    titlestring2=titlestring+".png";
    %saveas(figure_k,titlestring2);
    tempstring="print  -dpng "+titlestring2;
    eval(tempstring); 
elseif wt_mode ==2
    disp('for matlab 3D plot');
    num_k = sqrt(length(kplane));
    EIGENCAR = bandlist;
    [nband,~] = size(EIGENCAR);
    kx_mesh = reshape(kplane(:,1),num_k,num_k);
    ky_mesh = reshape(kplane(:,2),num_k,num_k);
    for Ei = 1:nband
        E_mesh = reshape(EIGENCAR(Ei,:),num_k,num_k);
        surf(kx_mesh,ky_mesh,E_mesh,'Parent',axes1,'FaceColor','interp',...
            'EdgeColor','none','FaceAlpha','0.6','DisplayName',num2str(Ei));
        hold on
    end
    if contourf_mode ==1
        for  Ei = 1:nband
        E_mesh = reshape(EIGENCAR(Ei,:),num_k,num_k);
        contourf(kx_mesh,ky_mesh, E_mesh,4)
        hold on

        end
    end
    xlabel('k_x ','Fontname','Helvetica')
    ylabel('k_y','Fontname','Helvetica')
    zlabel('Energy (eV)')
    axis image
    axis vis3d
    shading interp;
    colormap(hsv);
    set(gca,'FontSize',30);
    zlim(Ecut);
    view(axes1,[-72 17]);
    
end
    
end