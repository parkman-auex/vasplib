function c = redbluek(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [0 0 0], then [0 0 0] to [1 0 0];
    m1 = m*0.5;
    r = 0*ones(m1,1);
    g = 0*ones(m1,1);
    b = (m1-1:-1:0)'/max(m1-1,1);
    r = [r;(0:m1-1)'/max(m1-1,1)];
    g = [g;0*ones(m1,1)];
    b = [b;0*ones(m1,1)];
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = 0*ones(m1,1);
    g = 0*ones(m1,1);
    b = (m1-1:-1:0)'/max(m1-1,1);
    r = [r;(0:m1-1)'/max(m1-1,1)];
    g = [g;0*ones(m1,1)];
    b = [b;0*ones(m1,1)];
end

c = [r g b]; 
