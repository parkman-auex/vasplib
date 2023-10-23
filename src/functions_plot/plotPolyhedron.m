function varargout = plotPolyhedron(V,color,alpha,options)
arguments
    V = [];
    color = 'r';
    alpha = 0.2;
    options.ax = handle([]);
    options.OriginPoint = [0 0 0];
    options.Orientation = eye(3);
end

if isempty(options.ax)
    Fig =  Figs(1,1);
    ax = Fig.axes(1);
else
    ax = options.ax;
end
% if nargin < 2
%     color = 'r';
% end
% if nargin <3
%     alpha = 0.2;
% end
% if nargin < 4
%     [fig,ax] = creat_figure();
% end
if size(V,2) == 3
    %
%     V = V*options.Orientation;
    V = V + options.OriginPoint;
    %
    % https://zhuanlan.zhihu.com/p/357836869
    x = V(:,1);
    y = V(:,2);
    z = V(:,3);
    
    shp=alphaShape(x,y,z);
    shp.Alpha=20*max([max(x)-min(x),max(y)-min(y),max(z)-min(z)])/2;
    [elements,nodes]=boundaryFacets(shp);
    normalvec=zeros(size(elements));
    for k=1:size(elements,1)
        normalvec(k,:)=cross(nodes(elements(k,2),:)-nodes(elements(k,1),:),nodes(elements(k,3),:)-nodes(elements(k,1),:));
    end
    faceNear=nchoosek(1:size(elements,1),2);
    trisurf(elements,nodes(:,1),nodes(:,2),nodes(:,3),'FaceColor',color,'FaceAlpha',alpha,'EdgeColor','none','Parent',ax);
    for k=1:size(faceNear,1)
        isEdge=ismember(elements(faceNear(k,1),:),elements(faceNear(k,2),:));
        if sum(isEdge)>1 && vecnorm(cross(normalvec(faceNear(k,1),:),normalvec(faceNear(k,2),:)),1)/(vecnorm(normalvec(faceNear(k,1),:),1)*vecnorm(normalvec(faceNear(k,2),:),1))>1e-6
            plot3(ax,nodes(elements(faceNear(k,1),isEdge),1),nodes(elements(faceNear(k,1),isEdge),2),nodes(elements(faceNear(k,1),isEdge),3),'Color','k');
        end
    end
    axis(ax,'equal');
    view(ax,3);
elseif size(V,2) == 2
    if isequal(options.OriginPoint,[0 0 0]) && isequal(options.Orientation ,eye(3))
        % https://zhuanlan.zhihu.com/p/357836869
        x = V(:,1);
        y = V(:,2);
        pgon=polyshape(x,y);
        %shp.Alpha=20*max([max(x)-min(x),max(y)-min(y)])/2;
        %[elements,nodes]=boundaryFacets(shp);
        %    normalvec=zeros(size(elements));
        %     for k=1:size(elements,1)
        %         normalvec(k,:)=cross(nodes(elements(k,2),:)-nodes(elements(k,1),:),nodes(elements(k,3),:)-nodes(elements(k,1),:));
        %     end
        %faceNear=nchoosek(1:size(elements,1),2);
        %hold on
        %polyshape(nodes(:,1),nodes(:,2),'FaceColor',color,'FaceAlpha',alpha,'EdgeColor','none','Parent',ax);
        %
        plot(ax,pgon,'FaceColor',color,'FaceAlpha',alpha);
        axis(ax,'equal');
        view(ax,2);
    else
        V = V *options.Orientation;
        V = V + options.OriginPoint;
        patch(ax,V(:,1),V(:,2),V(:,3),'FaceColor',color,'FaceAlpha',alpha);
        axis(ax,'equal');
        view(ax,3);
    end

else
    
end

    %-------- return --------
    if nargout  == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    end
    if nargout  == 1
        varargout{1} = ax;
    end
end