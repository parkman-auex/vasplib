classdef gnuplot
    %COLORMAP from gnuplot
    properties(Constant)
        traditional = ColorMap.gnuplot.colormapformulae('traditional');
        rainbow = ColorMap.gnuplot.colormapformulae('rainbow');
        ocean = ColorMap.gnuplot.colormapformulae('ocean');
        pm3d = ColorMap.gnuplot.colormapformulae('pm3d'); 
        
    end
    
    methods
        function cmap = gnuplot()
            cmap = ColorMap.gnuplot.colormapformulae('traditional');
        end 
    end
    methods(Static)%get
    end
    methods(Static)
        function c = colormapformulae( type, m )
            if( ~exist('type','var') || isempty(type) )
                type = [7,5,15];
            end
            
            if( ~exist('m','var') || isempty(m) )
                m =256;
            end
            
            if( strcmpi(type, 'parula') )
                c = parula(m);
            elseif( strcmpi(type, 'jet') )
                c = jet(m);
            elseif( strcmpi(type, 'hsv') )
                c = hsv(m);
            elseif( strcmpi(type, 'cool') )
                c = cool(m);
            elseif( strcmpi(type, 'spring') )
                c = spring(m);
            elseif( strcmpi(type, 'autumn') )
                c = autumn(m);
            elseif( strcmpi(type, 'winter') )
                c = winter(m);
            elseif( strcmpi(type, 'bone') )
                c = bone(m);
            elseif( strcmpi(type, 'copper') )
                c = copper(m);
            elseif( strcmpi(type, 'pink') )
                c = pink(m);
                
            else
                if( isa(type, 'char') )
                    type = ColorMap.gnuplot.type2ids( type );
                end
                
                x = [0:1.0/(m-1):1];
                c = zeros(numel(x),3);
                for i=1:3
                    c(:,i) = ColorMap.gnuplot.formulae(x,type(i));
                end
                c(c<0)=0;
                c(c>1)=1;
            end
            
        end
        function y = formulae( x, id )
            % It is color transform function
            %
            % Inputs:
            %  x: range of input is [0,1]
            %  id: -36 to 36
            %
            % Outputs:
            %  y: output
            %
            % functions:
            %   0: 0
            %   1: 0.5
            %   2: 1
            %   3: x
            %   4: x^2
            %   5: x^3
            %   6: x^4
            %   7: sqrt(x)
            %   8: sqrt(sqrt(x))
            %   9: sin(90x)
            %  10: cos(90x)
            %  11: |x-0.5|
            %  12: (2x-1)^2
            %  13: sin(180x)
            %  14: |cos(180x)|
            %  15: sin(360x)
            %  16: cos(360x)
            %  17: |sin(360x)|
            %  18: |cos(360x)|
            %  19: |sin(720x)|
            %  20: |cos(720x)|
            %  21: 3x
            %  22: 3x-1
            %  23: 3x-2
            %  24: |3x-1|
            %  25: |3x-2|
            %  26: (3x-1)/2
            %  27: (3x-2)/2
            %  28: |(3x-1)/2|
            %  29: |(3x-2)/2|
            %  30: x/0.32-0.78125
            %  31: 2*x-0.84
            %  32: 4x;1;-2x+1.84;x/0.08-11.5
            %  33: |2*x - 0.5|
            %  34: 2*x
            %  35: 2*x - 0.5
            %  36: 2*x - 1
            % negative numbers mean inverted (1-output)
            %
            % version:201800402
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % colormap formulae:                                       %
            %                                                          %
            % Copyright (C) 2018 Masayuki Tanaka. All rights reserved. %
            %                    mtanaka@sc.e.titech.ac.jp             %
            %                                                          %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % http://juluribk.com/2011/05/18/843/
            % http://gnuplot-tricks.blogspot.jp/2009/06/comment-on-phonged-surfaces-in-gnuplot.html
            
            x(x<0)=0;
            x(x>1)=1;
            aid = abs(id);
            
            if( aid == 0 )
                y = zeros( size(x) );
                
            elseif( aid == 1 )
                y = ones( size(x) ) * 0.5;
                
            elseif( aid == 1 )
                y = ones( size(x) ) * 0.5;
                
            elseif( aid == 2 )
                y = ones( size(x) );
                
            elseif( aid == 3 )
                y = x;
                
            elseif( aid == 4 )
                y = x .^ 2;
                
            elseif( aid == 5 )
                y = x .^ 3;
                
            elseif( aid == 6 )
                y = x .^ 4;
                
            elseif( aid == 7 )
                y = sqrt(x);
                
            elseif( aid == 8 )
                y = x .^ 0.25;
                
            elseif( aid == 9 )
                y = sin( 0.5 * pi * x );
                
            elseif( aid == 10 )
                y = cos( 0.5 * pi * x );
                
            elseif( aid == 11 )
                y = abs(x-0.5);
                
            elseif( aid == 12 )
                y = (2*x-1) .* (2*x-1);
                
            elseif( aid == 13 )
                y = sin( pi * x );
                
            elseif( aid == 14 )
                y = abs( cos( pi * x ) );
                
            elseif( aid == 15 )
                y = sin( 2. * pi * x );
                
            elseif( aid == 16 )
                y = cos( 2. * pi * x );
                
            elseif( aid == 17 )
                y = abs( sin( 2. * pi * x ) );
                
            elseif( aid == 18 )
                y = abs( cos( 2. * pi * x ) );
                
            elseif( aid == 19 )
                y = abs( sin( 4. * pi * x ) );
                
            elseif( aid == 20 )
                y = abs( cos( 4. * pi * x ) );
                
            elseif( aid == 21 )
                y = 3 * x;
                
            elseif( aid == 22 )
                y = 3 * x - 1;
                
            elseif( aid == 23 )
                y = 3 * x - 2;
                
            elseif( aid == 24 )
                y = abs( 3 * x - 1);
                
            elseif( aid == 25 )
                y = abs( 3 * x - 2);
                
            elseif( aid == 26 )
                y = (3*x-1)/2.0;
                
            elseif( aid == 27 )
                y = (3*x-2)/2.0;
                
            elseif( aid == 28 )
                y = abs( (3*x-1)/2.0 );
                
            elseif( aid == 29 )
                y = abs( (3*x-2)/2.0 );
                
            elseif( aid == 30 )
                y = x/0.32-0.78125;
                
            elseif( aid == 31 )
                y = 2*x-0.84;
                
            elseif( aid == 32 )
                y = 4 * x;
                y(x>0.25 & x<=0.5) = 1;
                y(x>0.5 & x<=0.75) = -2 * x(x>0.5 & x<=0.75) + 1.84;
                y(x>0.75) = x(x>0.75) / 0.08 - 11.5;
                
            elseif( aid == 33 )
                y = abs( 2*x - 0.5 );
                
            elseif( aid == 34 )
                y = 2*x;
                
            elseif( aid == 35 )
                y = 2*x-0.5;
            elseif( aid == 36 )
                y = 2*x-1;
            else
                error('Invalidate input for id. id should be -36 to 36.');
            end
            
            if( id < 0 )
                y = 1-y;
            end
            
        end
        function col = gray2pseudocolor( x, type )
            % It converts gray image to pseudocolor image
            %
            % Inputs:
            %  x: gray image whose range is [0,1]
            %  type: colormap type
            %
            % Output:
            %  col: pseudo color image [0,1]
            %
            % typical colormap names of pnm3d:
            %  [ 7, 5,15]:traditional,pm3d,black-blue-red-yellow
            %  [ 3,11, 6]:green-red-violet
            %  [23,28, 3]:ocean,green-blue-white
            %  [21,22,23]:hot,black-red-yellow-white
            %  [30,31,32]:color printable on gray,black-blue-violet-yellow-white
            %  [33,13,10]:rainbow,blue-green-yellow-red
            %  [34,35,36]:AFM hot,black-red-yellow-white
            %
            % functions:
            %   0: 0
            %   1: 0.5
            %   2: 1
            %   3: x
            %   4: x^2
            %   5: x^3
            %   6: x^4
            %   7: sqrt(x)
            %   8: sqrt(sqrt(x))
            %   9: sin(90x)
            %  10: cos(90x)
            %  11: |x-0.5|
            %  12: (2x-1)^2
            %  13: sin(180x)
            %  14: |cos(180x)|
            %  15: sin(360x)
            %  16: cos(360x)
            %  17: |sin(360x)|
            %  18: |cos(360x)|
            %  19: |sin(720x)|
            %  20: |cos(720x)|
            %  21: 3x
            %  22: 3x-1
            %  23: 3x-2
            %  24: |3x-1|
            %  25: |3x-2|
            %  26: (3x-1)/2
            %  27: (3x-2)/2
            %  28: |(3x-1)/2|
            %  29: |(3x-2)/2|
            %  30: x/0.32-0.78125
            %  31: 2*x-0.84
            %  32: 4x;1;-2x+1.84;x/0.08-11.5
            %  33: |2*x - 0.5|
            %  34: 2*x
            %  35: 2*x - 0.5
            %  36: 2*x - 1
            % negative numbers mean inverted (1-output)
            %
            % version:201800402
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % colormap formulae:                                       %
            %                                                          %
            % Copyright (C) 2018 Masayuki Tanaka. All rights reserved. %
            %                    mtanaka@sc.e.titech.ac.jp             %
            %                                                          %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if( ~exist('type','var') || isempty(type) )
                type = [7,5,15];
            end
            
            if( isa(type, 'char') )
                type = ColorMap.gnuplot.type2ids( type );
            end
            
            col = zeros( size(x,1), size(x,2), 3 );
            
            for i=1:3
                col(:,:,i) = ColorMap.gnuplot.formulae( x, type(i) );
            end
            
            col(col<0) = 0;
            col(col>1) = 1;
            
        end
        function ids = type2ids(type )
            % It converts string to ids
            %
            % Input:
            %  type: name
            %
            % Output:
            %  ids: three dimensional id
            %
            % typical colormap names of pnm3d:
            %  [ 7, 5,15]:traditional,pm3d,black-blue-red-yellow
            %  [ 3,11, 6]:green-red-violet
            %  [23,28, 3]:ocean,green-blue-white
            %  [21,22,23]:hot,black-red-yellow-white
            %  [30,31,32]:color printable on gray,black-blue-violet-yellow-white
            %  [33,13,10]:rainbow,blue-green-yellow-red
            %  [34,35,36]:AFM hot,black-red-yellow-white
            %  [ 3, 3, 3]:gray
            %
            % version:201800402
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % colormap formulae:                                       %
            %                                                          %
            % Copyright (C) 2018 Masayuki Tanaka. All rights reserved. %
            %                    mtanaka@sc.e.titech.ac.jp             %
            %                                                          %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if( ~exist('type','var') || isempty(type) )
                fprintf('[ 7, 5,15]:traditional,pm3d,black-blue-red-yellow\n');
                fprintf('[ 3,11, 6]:green-red-violet\n');
                fprintf('[23,28, 3]:ocean,green-blue-white\n');
                fprintf('[21,22,23]:hot,black-red-yellow-white\n');
                fprintf('[30,31,32]:color printable on gray,black-blue-violet-yellow-white\n');
                fprintf('[33,13,10]:rainbow,blue-green-yellow-red\n');
                fprintf('[34,35,36]:AFM hot,black-red-yellow-white\n');
                fprintf('[ 3, 3, 3]:gray\n');
                
                
            else
                
                if( strcmpi( type, 'traditional' ) || strcmpi( type, 'pm3d' ) || strcmpi( type, 'black-blue-red-yellow' ) )
                    ids = [7,5,15];
                elseif( strcmpi( type, 'green-red-violet' ) )
                    ids = [3,11,6];
                elseif( strcmpi( type, 'ocean' ) || strcmpi( type, 'green-blue-white' ) )
                    ids = [23,28,3];
                elseif( strcmpi( type, 'hot' ) || strcmpi( type, 'black-red-yellow-white' ) )
                    ids = [21,22,23];
                elseif( strcmpi(type, 'color printable on gray') || strcmpi( type, 'black-blue-violet-yellow-white' ) )
                    ids = [30,31,32];
                elseif( strcmpi(type, 'rainbow') || strcmpi( type, 'blue-green-yellow-red' ) )
                    ids = [33,13,10];
                elseif( strcmpi(type, 'AFM hot') || strcmpi( type, 'black-red-yellow-white' ) )
                    ids = [34,35,36];
                elseif( strcmpi(type, 'gray') )
                    ids = [3,3,3];
                else
                    error( 'Unknown type' );
                    % set palette defined ( -10 "#194eff", 0 "white", 10 "red" )
                end
                
            end
        end
    end
end

