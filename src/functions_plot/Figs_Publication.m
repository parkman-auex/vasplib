function [fig, axes] = Figs_Publication(m,n)
arguments
    m = 1
    n = 1
end
%%
ax_size = [7 7];
ax_bias = [2.8 1.8];
fig_size = [10.5 9.9]; % 3x2 A4 size
fig_bias = [5 5];
%%
fig = figure();
set(fig,'unit', 'centimeters', 'position', [fig_bias fig_size(1)*n fig_size(2)*m], 'Color', 'w');

p = 0;
% alphabet_list = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"];

axes = gobjects(m,n);
for i = 1:m
    for j = 1:n
        p = p+1;
        axes(i,j) = subplot(m,n,p);
        set(axes(i,j), 'unit', 'centimeters', 'Position',...
            [ax_bias(1) + fig_size(1)*(j-1),...
             ax_bias(2) + fig_size(2)*(m-i),...
             ax_size]);
        set(axes(i,j), 'FontSize', 18, 'FontName', 'Arial');
        set(axes(i,j), 'LineWidth',1,'Box','on')
        % title(axes(i,j), "("+alphabet_list(p)+")", 'Position',[-0.32, 1.05], 'Units','centimeters')
    end
end
end