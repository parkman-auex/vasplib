% plotg creates line plots with color gradients 
%
% syntax
%   plotg(y)
%   plotg(x,y)
%   plotg(x,y,z)
%   plotg(x,y,z,cdata)
%   plotg('Property', 'Value', ...)
%   h = plotg(_)
%
% Description
%   plotg(y) plots the columns of y versus their index with color gradient
%   based on 'parula' colormap that varies from the first to last index of 
%   y vector.
%
%   plot(x,y) plots vector y versus vector x with 'parula' colormap that 
%   varies from the first to last indices of x-y pairs.
%
%   plot(x,y,z) plots a line in 3-space through the points whose 
%   coordinates are indicated by x, y, z vectors of the same size. The
%   color of the line changes from the first to the last coordinate based
%   on 'parula' colormap
%
%   plot(x,y,z, cdata) plots a line in 3-space through the points whose 
%   coordinates are indicated by x, y, z vectors of the same size. The
%   color at each point of the line is defined based on cdata as
%
%       if cdata is a colormap name (or a name of function that creates
%       colormaps) as a string, then plotg will use the corresponding 
%       colormap to define the colors of the line from its first index to 
%       the last.
%
%       if cdata is a n by 3 vector with values between 0 and 1, then plotg
%       uses cdata as a matrix of rgb values. The color at each point of
%       the line is defined by interpolating the values in cdata such that
%       the color of the first and last points on the line match with 
%       cdata(1,:) and cdata(end,:), respectively.
%
%       if cdata is a vector with the same length as x,y and/or z vectors,
%       then plotg will assign as the values of the points along the line
%       and the figure colormap is used to defined the color of the points
%       (based on their values).
%
%   h = plotg(_) returns the handle to the plotg object. The properties of
%   the plotg object are
%       Handle: the handle to patch used to generate plotg
%       XData
%       YData
%       ZData
%       CData
%       LineWidth
%       LineStyle
%       
%
% Examples
%       >> t = linspace(0,2*pi,100);
%       >> x = cos(t);
%       >> y = sin(t);
%       >> z = t;
%       >> h = plotg(x,y,z,'jet');
%       >> view(3);
%
%  Similar result can be obtained by defining the properties of plotg
%       >> h = plotg('XData',x,'YData',y,'ZData',z,'CData','jet');
%       >> view(3);
%
%  The properties of plotg can also be modified afterwards
%       >> set(h,'LineWidth',3,'LineStyle','--');
%
%  For other examples please check plotgExample.m
%__________________________________________________________________________
% Copyright and disclaimer
% This software (plotg.m) is provided by the provider (Siamak G. Faal) 
% "as is" and “with all faults.” The provider makes no representations or 
% warranties of any kind concerning the safety, suitability, lack of 
% viruses, inaccuracies, typographical errors, or other harmful components 
% of this software product. There are inherent dangers in the use of any 
% software, and you are solely responsible for determining whether this 
% software product is compatible with your application, equipment and other
% software installed on your equipment. You are also solely responsible for
% the protection of your equipment and backup of your data, and the 
% provider will not be liable for any damages you may suffer in connection 
% with using, modifying, or distributing this software product. 
% Additionally, all data, information and examples provided in this 
% document are for informational purposes only. 
% All information is provided on an “as is” basis and the author 
% (Siamak G. Faal) makes no representations as to accuracy, completeness, 
% correctness, suitability, or validity of any information provided. 
% The author will not be liable for any errors, omissions, or delays in 
% this information or any losses, injuries, or damages arising from its 
% use. Your use of any information and/or examples is entirely at your own 
% risk. Should the software/information/examples prove defective, you 
% assume the entire cost of all service, repair or correction.
%
% Author: Siamak G. Faal
%         sghorbanifaal@wpi.edu
%         https://users.wpi.edu/~sghorbanifaal/
%
% Version: 1.5
% Latest update: September 16, 2017
%__________________________________________________________________________
classdef plotg <  handle & matlab.mixin.SetGet
    properties (SetObservable)
        Handle = [];
        XData = [];
        YData = [];
        ZData = [];
        CData = 'parula';
        LineWidth = 1;
        LineStyle = '-';
    end
    
    %______________________________________________________________________
    % Methods - Private
    methods(Access = private)
        %Attach listeners to the visualObject properties
        function attachListener(obj)
            addlistener(obj,'XData','PostSet',@(metaProp,eventData)plotg.dataUpdate(metaProp,eventData,obj));
            addlistener(obj,'YData','PostSet',@(metaProp,eventData)plotg.dataUpdate(metaProp,eventData,obj));
            addlistener(obj,'ZData','PostSet',@(metaProp,eventData)plotg.dataUpdate(metaProp,eventData,obj));
            addlistener(obj,'CData','PostSet',@(metaProp,eventData)plotg.dataUpdate(metaProp,eventData,obj));
            addlistener(obj,'LineWidth','PostSet',@(metaProp,eventData)plotg.appearanceChange(metaProp,eventData,obj));
            addlistener(obj,'LineStyle','PostSet',@(metaProp,eventData)plotg.appearanceChange(metaProp,eventData,obj));
        end
    end
        
    methods    
        function obj = plotg(varargin)
            % obj = plotg(YData)
            % obj = plotg(XData, YData)
            % obj = plotg(XData, YData, ZData)
            % obj = plotg(XData, YData, ZData, CData)
            % obj = plotg('Property', 'Value', ...)
            obj.attachListener
            switch nargin
                case 0
                case 1 
                    obj.XData = 1:length(varargin{1});
                    obj.YData = varargin{1};
                case 2
                    obj.XData = varargin{1};
                    obj.YData = varargin{2};
                case 3
                    obj.XData = varargin{1};
                    obj.YData = varargin{2};
                    obj.ZData = varargin{3};
                case 4
                    obj.XData = varargin{1};
                    obj.YData = varargin{2};
                    obj.ZData = varargin{3};
                    obj.CData = varargin{4};
                otherwise
                    for n = 1:2:nargin-1
                        obj.(varargin{n}) = varargin{n+1};
                    end
            end
        end
    end
    
        %______________________________________________________________________
    % Methods - Static
    methods(Static)
        %__________________________________________________________________
        % shapeChange - A listener for: XData, YData, ZData, CData
        function dataUpdate(~,~,obj)            
            if(isempty(obj.Handle))
                obj.Handle = patch('Faces',[],'Vertices',[],...
                    'FaceColor','none','EdgeColor','interp');
            end
            
            n = length(obj.XData);
            F = [1:n n-1:-1:1];
            V = get(obj.Handle,'Vertices');
            colors = creatColorMap(obj.CData,n);
            if(length(obj.XData)==length(obj.YData))
                if(isempty(obj.ZData))
                    V = [reshape(obj.XData,n,1) reshape(obj.YData,n,1)];
                elseif(length(obj.XData)==length(obj.ZData))
                    V = [reshape(obj.XData,n,1) reshape(obj.YData,n,1) reshape(obj.ZData,n,1)];
                end
                set(obj.Handle,'Vertices',V,'Faces',F,'FaceVertexCData',colors);
            end
        end
        %__________________________________________________________________
        % appearanceChange - A listener for: LineWidth, LineStyle
        function appearanceChange(metaProp,eventData,obj) 
            h = eventData.AffectedObject;
            propName = metaProp.Name;
            set(obj.Handle,propName,h.(propName));
        end 
        function about() 
            fprintf(' plotg.m \n Version 1.5 \n by Siamak Faal\n');
        end
    end
end

function colors = creatColorMap(CData,n)

if(ischar(CData))
    CData = str2func(CData);
    f = @(n)CData(n);
    colors = f(n);
elseif(size(CData,2) == 3)
    colorMax = max(max(CData));
    colorMin = min(min(CData));

    if(colorMax > 1 || colorMin < 0)
        error('Colormap must have values in [0,1].');
    end

    if(n==1)
        colors = CData(1,:);
    else

        x = linspace(1,n,size(CData,1));
        p = 1:n;

        R = spline(x,CData(:,1),p);
        G = spline(x,CData(:,2),p);
        B = spline(x,CData(:,3),p);

        l = min([min(R), min(G), min(B)]);
        h = max([max(R), max(G), max(B)]);

        colors = ([reshape(R,n,1), reshape(G,n,1), reshape(B,n,1)] - l)/(h-l);
    end
    
elseif(length(CData) == n)
    colorMap = get(gcf,'Colormap');
    nC = size(colorMap,1)-1;
    data = (CData - min(CData));
    data = nC*data/max(data)+1;
    colors = zeros(n,3);
    for i=1:n
        l = floor(data(i));
        u = ceil(data(i));
        colors(i,:) = (colorMap(u,:)-colorMap(l,:))/(u-l)*(data(i)-l)+...
            colorMap(l,:);
    end
else
    error('Invalid color specification input.');
end
end