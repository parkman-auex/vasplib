classdef Figs < matlab.graphics.chartcontainer.ChartContainer
    properties
        rows double = 1
        cols double = 1
        axes(:,:)
        fig_handle
        
        TitleText string = ""
        XLabelText string = ""
        YLabelText string = ""
        XLim double = []
        YLim double = []
        XTicks
        XTickLabels
        
        FontName = 'Helvetica'
        FontSize = 24   
        LineWidth = 1
        Color = [0.9400 0.9400 0.9400]
    end

    methods
        function obj = Figs(rows_in,cols_in)
            % fig = figure();
            args = {'rows',rows_in,'cols',cols_in};
            obj = obj@matlab.graphics.chartcontainer.ChartContainer(args{:});
           
            drawnow
        end              
    end
    
    methods(Access=protected)    
        function setup(~)
        end
        
        function update(obj)
            tcl = getLayout(obj);
            delete(tcl.Children);              
            tcl.GridSize = [obj.rows obj.cols];
            for i = 1:obj.rows * obj.cols
                nexttile(tcl,i);
            end

            % Populate the layout with the axes
            obj.axes = gobjects(obj.rows, obj.cols);
            for m = 1:obj.rows
                for n = 1:obj.cols
                    % Get the axes at the current row/column
                    t = (m-1) * obj.cols + n;
                    obj.axes(m,n) = nexttile(tcl,t);
                end
            end

         
            

            
            unit_base = [0.4,0.4];
            unit_factor = sqrt([obj.cols, obj.rows]);
            unit_base = (unit_base .* unit_factor);
            unit_base =  unit_base/max(unit_base)*0.6;
            set(gcf,'Units','normalized');            
            set(gcf,'Position',[0.1 0.1 unit_base(1) unit_base(2)]);
            
    
%             
%             % Chart style
%             % some properties are in super class and can be use directly.
% %             tcl.TileSpacing = 'compact';
% %             tcl.Padding = 'tight';            
% %             title(tcl,obj.TitleText,'FontSize',obj.FontSize,'FontName',obj.FontName);
%             
%             xlabel(obj.axes,obj.XLabelText);
%             ylabel(obj.axes,obj.YLabelText);
%             if ~isempty(obj.XLim)
%                 xlim(obj.axes,obj.XLim);
%             end
%             if ~isempty(obj.YLim)
%                 ylim(obj.axes,obj.YLim);
%             end
%             xticks(obj.axes,obj.XTicks);
%             xticklabels(obj.axes,obj.XTickLabels);
%             set(obj.axes,'FontSize',obj.FontSize,'FontName',obj.FontName,...
%                 'LineWidth',obj.LineWidth);
            % hold(obj.axes,'on');
            % box(obj.axes,'on');
            %set(obj.Parent,'Color','w');
        end
    end
end

