function [Hopping_handle] = Plot_Hopping(vectorL,Hoppinglist,Rm,orbL,options)
    arguments
        vectorL;
        Hoppinglist;
        Rm;
        orbL;
        options.ax = gca();
        options.cmap = @turbo;
        options.TwoD = false;
        options.W = 0.5;
        options.H = 1.5;
        options.LineWidth = 1;
        options.norm = false
    end
    
    if options.TwoD
        Rm(3,:) = [0,0,1];
    end
    ax = options.ax ;
    vectorList = double(vectorL);
    if options.norm
        Hoppinglist = abs(Hoppinglist);
    end
    [UniqueList,~,ic] = unique(Hoppinglist,'rows');
    ColorList = options.cmap(length(UniqueList));
    for i = 1:size(vectorList,1)
        % $H_{i j}^{\mathbf{k}}=\left\langle\chi_{i}^{\mathbf{k}}|H| \chi_{j}^{\mathbf{k}}\right\rangle=\sum_{\mathbf{R}} e^{i \mathbf{k} \cdot\left(\mathbf{R}+\mathbf{t}_{j}-\mathbf{t}_{i}\right)} H_{i j}(\mathbf{R})$
        Start = (orbL((vectorList(i,5)),:)+vectorList(i,1:3))*Rm;
        End = (orbL((vectorList(i,4)),:))*Rm;
        DisplayName = mat2str(vectorList(i,:));
        if isequal(Start,End)
            continue;
        end
        Arrow = Plot_arrow(Start,End...
            ,'ax',ax...
            ,'Color',ColorList(ic(i),:) ...
            ,'LineWidth',options.LineWidth ...
            ,'W',options.W ...
            ,'H',options.H ...
            ,'DisplayName',DisplayName...
            );
        Hopping_handle(i,:) = Arrow.';
    end
    axis(ax,'equal');
end
