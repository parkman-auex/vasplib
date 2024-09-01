function c = redbluelight(m)
arguments
    m = 101
end

c= ColorMap.customcolormap(linspace(0,1,7),...
    {'#aa2409','#f64623','#fea290','#ffffff','#92c7fc','#258cf4','#0d59a5'},m);

end