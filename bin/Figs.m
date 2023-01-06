function [f_struct,Ax] = Figs(rows, cols,opts)
% return a struct contains the handle and the axes of a tiledlayout
% f_struct.fig = handle
% f_struct.axes = a rows * cols GraphicsPlaceholder array
arguments
    rows double = 1
    cols double = 1
    opts.FontName = 'Helvetica';
    opts.FontSize = 24;
end
%%
figure();
tcl = tiledlayout(rows,cols);
tcl.TileSpacing = 'compact';
tcl.Padding = 'tight';            
f_struct.fig = gcf;
for i = 1 : rows*cols
    nexttile(tcl,i);
end

% Populate the layout with the axes
f_struct.axes = gobjects(rows, cols);
for m = 1:rows
    for n = 1:cols
        % Get the axes at the current row/column
        t = (m-1) * cols + n;
        f_struct.axes(m,n) = nexttile(tcl,t);
    end
end
%%
unit_base = [0.4,0.4];
unit_factor = sqrt([cols, rows]);
unit_base = (unit_base .* unit_factor);
unit_base =  unit_base/max(unit_base)*0.6;
set(f_struct.fig,'Units','normalized');            
set(f_struct.fig,'Position',[0.1 0.1 unit_base(1) unit_base(2)]);         
%%
FontName = opts.FontName;
FontSize = opts.FontSize;
LineWidth = 1;
% Color = [0.9400 0.9400 0.9400];
set( f_struct.axes,'FontSize',FontSize,'FontName',FontName,'LineWidth',LineWidth);
box( f_struct.axes,'on');
set( f_struct.fig,'Color','w');
hold(f_struct.axes,'on');
Ax = f_struct.axes;
end