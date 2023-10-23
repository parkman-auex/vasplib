function Arrow = Plot_arrow(Start,End,options)
arguments
    Start = [0,0,0];
    End = [1,1,1];
    options.ax = gca();
    options.Color = 'b';
    options.LineStyle  = '*';
    options.LineWidth = 1.0;
    options.W = 1;
    options.H = 2;
    options.DisplayName = 'Hopping';
    options.CData = [];
end
    import vasplib_tool_outer.*;
    S = [ 'd',options.LineStyle , char(string(options.LineWidth))];
    W = options.W;
    H = options.H;
    Arrow = vasplib_tool_outer.arrow3(Start,End,S,W,H);
    %     Arrow =  vasplib_tool_outer.arrow3('update');
    Line = handle(Arrow(1));
    Surf = handle(Arrow(2));
    Line.DisplayName = options.DisplayName;
    Line.Color = options.Color;
    Surf.DisplayName = options.DisplayName;
    Surf.CData(:,:,1) = Line.Color(1);
    Surf.CData(:,:,2) = Line.Color(2);
    Surf.CData(:,:,3) = Line.Color(3);
    if ~isempty(options.CData)
        Surf.CData = options.CData;
    end
end