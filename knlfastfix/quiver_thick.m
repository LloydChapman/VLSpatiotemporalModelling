function varargout = quiver_thick(x1,y1,u,v,varargin)
%quiver_thick plot fancy arrows like arrow and quiver
%__________________________________________________________________________
%Introduction
%
%QUIVER_THICK allows users to plot arrows just like the syntax of the
%functions 'arrow' and 'quiver' and with identical speeds (so practically
%instantaneous), though as fancy as (slower) curvec plots (straight arrows
%though). It also includes functionality to plot magnitudes using colors in
%the filled arrows (lacking functionality in curvec) and supports
%(ir)regular grids, meshes and simple vectors and has a lot of user
%flexibility by using optional <keyword,value> pairs.
%
%Have a look at the help below to get an idea of the syntax & functionality
%__________________________________________________________________________
%Help
%
%quiver_thick uses start locations (X,Y) and components (U,V) to generate 2D
%plots of arrows as filled patches. Thus, the function works just like the
%'quiver' function. The locations and components (input) can be specified 
%on an (irregular) grid/mesh or basic matlab vectors. When using start (x1,y1)
%and end (x2,y2) locations, one is referred to arrow_thick, which works
%just like the well known function 'arrow' with (x1,y1,x2,y2). By default,
%arrows are colored (filled) according to their magnitude.
%
%quiver_thick is actually a wrapper for arrow_thick, it therefor uses identical
%additional keywords and value pairs. Though, three additional keyword-value
%pairs have been added ('scaling', 'scaling_type' & 'uniform_length') that
%are relevant for quiver type plotting (using vectorial components)
%
%Syntax:
%
%<output> = quiver_thick(X,Y,U,V,<keyword,value>);
%
%Apart from X, Y, U & V please note that all other in- and output
%variables (depicted by <...>) are optional
%
%Input variables:
%
%X                         Required: Start X-coordinates of arrows, can be
%                          any vector or matrix [MxN] M,N = 1:any integer
%Y                         Required: Start Y-coordinates of arrows, can be
%                          any vector or matrix [MxN] M,N = 1:any integer
%U                         Required: U-components (X-direction), can be
%                          any vector or matrix [MxN] M,N = 1:any integer
%V                         Required: V-components (Y-direction), can be
%                          any vector or matrix [MxN] M,N = 1:any integer
%
%All <keyword,value> pairs below are optional and shown as <name,default>
%
%Simply call quiver_thick to get a list of optional keywords and their defaults
%
%<'scaling',1>  
%                          'scaling' is the factor by which the length of
%                          the arrow is scaled. By default, an arrow with a
%                          U-component of 10 and a V-component of 0 will
%                          have a length of 10. When setting e.g. 'scaling'
%                          to 10, the length will be plotted as 100 (the
%                          components are scaled by a factor 'scaling'
%
%<'scaling_type','magnitude'>  
%                          By setting the 'scaling_type' to 'uniform' all
%                          plotted arrows have identical lengths, which
%                          means that only location and direction are
%                          shown. By default, the length of arrows depends
%                          on the magnitude of the U- and V-components and
%                          'scaling_type' is therefor 'magnitude' by default
%
%<'uniform_length',1>  
%                          When using 'uniform' for 'scaling_type' the
%                          uniform length of each arrow can be specified.
%                          By default, when using the 'uniform' option for
%                          'scaling_type', the length of each vector is 1.
%                          This can be changed by specifying the length of
%                          'uniform_length' as a number. This is for
%                          instance useful when plotting on geographic
%                          coordinates (by default an arrow would be 1 deg.,
%                          which is approx. 110 km along the equator)
%                          
%<'arrow_thickness',0.05>  
%                          The thickness of the arrow relative to its
%                          length. Can also be specified as an absolute
%                          value (by using the keyword below), the absolute
%                          value is based on the axis it is plotted in.
%                          When specifying a value of 0 (or lower) the
%                          arrow will be plotted as a line, just like the
%                          functions 'arrow' and 'quiver'. Default = 0.05
%
%<'arrow_thickness_option','relative'>  
%                          Switch for interpretation of the keyword
%                          arrow_thickness. Can be 'relative' (default) or
%                          'absolute' (as described above)'. Any other
%                          statement will trigger the default approach
%
%<'arrowhead_length',0.15>  
%                          The length of the arrowhead relative to its
%                          length. Can also be specified as an absolute
%                          value (by using the keyword below), the absolute
%                          value is based on the axis it is plotted in.
%                          When specifying a value of 0 (or lower) no
%                          arrowhead will be plotted, if any of the plotted
%                          arrowheads is longer than the arrow itself, the
%                          arrowhead is automatically limited to 90% of the
%                          arrows' length. Default = 0.15
%
%<'arrowhead_length_option','relative'>  
%                          Switch for interpretation of the keyword
%                          arrowhead_length. Can be 'relative' (default) or
%                          'absolute' (as described above)'. Any other
%                          statement will trigger the default approach
%
%<'max_arrowhead_angle',60>  
%                          Specification of the maximum arrowhead angle in
%                          in degrees, can be anywhere from 0 to 170 degrees
%
%<'plot_color','magnitude'>
%                          'plot_color' allows you to change the color of
%                          the area of the filled arrows area, by default,
%                          'magnitude' triggers different colors for all
%                          arrows based on their lengths, using the
%                          colormap as specified by the keyword 'colormap'.
%                          Other options are identical to those used in the
%                          function 'plot', e.g. 'r', 'g', 'b', 'c', 'm',
%                          'y', 'k', 'w' or any [r g b] vector (e.g. green
%                          [0 1 0]), which specifies a single color for all
%                          filled arrows
%
%<'plot_colorbar',1>
%                          When using 'plot_color' with 'magnitude'),
%                          a colorbar is automatically added to the figure.
%                          If you want to avoid this behaviour, simply set
%                          'plot_colorbar' to 0
%
%<'colormap',flipud(hot(64))>  
%                          When using 'plot_color' with 'magnitude' the
%                          colormap specified by 'colormap' is used. This
%                          can be any [Mx3] matrix (e.g. generated by jet,
%                          cool or any manually generated colormap). By
%                          default flipud(hot(64)) is used, which makes
%                          arrows seem 'hotter' when they are larger
%
%<'plot_to_current_figure',1>  
%                          When setting 'plot_to_current_figure' to 0, the
%                          arrows are not plotted in the current figure (or
%                          a new one if not present yet), when set to 1
%                          (which is default) the arrows are plotted
%
%<'axis_equal',1>  
%                          When setting 'axis_equal' to 0, the axis equal
%                          call within the function is skipped. This is
%                          only relevant if 'plot_to_current_figure' is 1.
%                          Do note that arrows can get distorted when not
%                          using axis_equal, as arrow settings are
%                          reflected upon the axis they are plotted on
%
%<'EdgeColor',[0 0 0]>
%                          The edges of the arrows are colored according to
%                          the color specified by 'EdgeColor', default ([0
%                          0 0]) means black, can also be specified by e.g.
%                          'r', 'g', 'k' apart from rgb values: [R G B]. To
%                          exclude the line around the arrows, set
%                          'LineStyle' to 'none'
%
%<'FaceAlpha',1>
%                          Defines the transparency of the arrows, can be
%                          any value from 0 (completely transparent) to 1
%                          (not transparent). Default value is 1
%                          
%<'LineWidth',0.5>
%                          Defines the width of the edge line, can be any
%                          positive value, default is 0.5, using 0 is not
%                          allowed, to exclude the line around the arrows,
%                          simply set 'LineStyle' to 'none'
%           
%<'LineStyle','-'>
%                          Defines the style of the line around the arrow,
%                          default is '-' (solid line), but can also be
%                          'none', '-.' or '--'
%
%
%Output variables (all ouput variables are optional):
%
%0 output variables        The ans will consist of the x and y values of
%                          the arrows in a structure, when plotted, a
%                          handle to the plotted arrows is also added to
%                          the structure
%1 output variable         The single output variable will consist of the x
%                          and y values of the arrows in a structure, when
%                          plotted, a handle to the plotted arrows is also
%                          added to the structure
%2 output variables        The two output variables will contain the x and
%                          y locations of the generated arrows in the two
%                          separate variables. Please note that a handle to
%                          the plotted arrows (if used) is not provided
%3 output variables        The first two output variables will contain the 
%                          x and y locations of the generated arrows in the 
%                          two separate variables. Please note that the
%                          handle to the plotted arrows (if used) is
%                          provided in the third argument. When a thrid
%                          argument is provided in case of plotting on, the
%                          thrid argument will be empty
%
%Note that in the output x- and y-locations the 2nd and 6th row contain the
%start and end locations or the arrows (also usable for arrow_thick)
%__________________________________________________________________________
%Examples
%_________
%Example 1 - All default settings:
%
%quiver_thick(zeros(1,36),...
%            zeros(1,36),...
%            (cosd(1:10:360)).*(36:-1:1),...
%            (sind(1:10:360)).*(36:-1:1));
%box on; axis tight; set(gca,'XTick',[],'YTick',[]);
%_________
%Example 2 - Making normal arrow/quiver plots, using both 'magnitude' and
%            'uniform' as arrow lengths ('scaling_type') in different plot boxes
%
%x=0:30; y=0:10; [X1,Y1]=meshgrid(x,y); X2=rand(11,31); Y2=2*rand(11,31);
%figure; subplot(2,1,1);
%[X_data_mag,Y_data_mag] = quiver_thick(X1,Y1,X2,Y2,...
%                 'arrow_thickness',0,...
%                 'arrowhead_length_option','absolute',...
%                 'scaling_type','magnitude',...
%                 'plot_color','k');
%box on; axis tight; set(gca,'XTick',[],'YTick',[]); title('Magnitude length');
%subplot(2,1,2);
%[X_data_uni,Y_data_uni] = quiver_thick(X1,Y1,X2,Y2,...
%                 'arrow_thickness',0,...
%                 'arrowhead_length_option','absolute',...
%                 'scaling_type','uniform',...
%                 'uniform_length',1,...
%                 'plot_color','k');
%box on; axis tight; set(gca,'XTick',[],'YTick',[]); title('Uniform lengths');
%_________
%Example 3 - Getting the arrow data and plot afterwards
%
%figure; plot(0:10,0:10,'k--'); hold on; plot(0:10,10:-1:0,'k--');
%OutputData = quiver_thick([6 6 4 4],[6 4 4 6],[3 3 -3 -3],[3 -3 -3 3],...
%                         'plot_to_current_figure',0);
%p1 = patch(OutputData.arrows_x_values,OutputData.arrows_y_values,'r');
%box on; axis equal; axis tight;
%set(p1,'FaceAlpha',0.6); set(gca,'XTick',[],'YTick',[]);
%_________
%Example 4 - Using many options
%
%X1=rand(999,1)*20-10; Y1=rand(999,1)*20-10; X2=0.2*X1; Y2=0.2*Y1;
%text(0,0,'KABOOM!','horizontalalignment','center',...
%         'fontsize',40,'fontweight','bold','color','r'); hold on;
%OutputData = quiver_thick(X1,Y1,X2,Y2,...
%            'arrow_thickness',0.1,...
%            'arrow_thickness_option','absolute',...
%            'colormap',hot(100),...
%            'FaceAlpha',0.8,...
%            'max_arrowhead_angle',40,...
%            'arrowhead_length',1,...
%            'arrowhead_length_option','absolute',...
%            'plot_colorbar',0);
%box on; axis tight; set(gca,'XTick',[],'YTick',[]);
%__________________________________________________________________________
%
%Contact Freek Scheel (freek.scheel@deltares.nl) if bugs are encountered
%              
%See also: arrow_thick patch arrow quiver curvec

%   --------------------------------------------------------------------
%   Copyright (C) 2014 Deltares
%       Freek Scheel
%       +31(0)88 335 8241
%       <freek.scheel@deltares.nl>;
%
%       Please contact me if errors occur.
%
%       Deltares
%       P.O. Box 177
%       2600 MH Delft
%       The Netherlands
%
%   This library is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this library.  If not, see <http://www.gnu.org/licenses/>.
%   --------------------------------------------------------------------

% Keywords from arrow_thick:
keywords.arrow_thickness         = 0.05;
keywords.arrow_thickness_option  = 'relative';

keywords.arrowhead_length        = 0.15;
keywords.arrowhead_length_option = 'relative';

keywords.max_arrowhead_angle     = 60;

keywords.plot_color              = 'magnitude';
keywords.colormap                = flipud(hot(64));
keywords.plot_colorbar           = 1;

keywords.EdgeColor               = [0 0 0];
keywords.EdgeAlpha               = 1;
keywords.FaceAlpha               = 1;
keywords.LineWidth               = 0.5;
keywords.LineStyle               = '-';

keywords.plot_to_current_figure  = 1;
keywords.axis_equal              = 1;

% Keywords from quiver_thick
keywords.scaling                 = 1;
keywords.scaling_type            = 'magnitude';
keywords.uniform_length          = 1;

if nargin == 0
    default_keywords = keywords
    return
end

if odd(length(varargin)) % Check for input consistency
    error('Odd number of input parameters for <keyword,value> pairs supplied, please check this');
end

% Set user defined changes to default keywords:

keywords = setproperty(keywords,varargin);

% Do some basic checks for input, all <keyword,value> pairs are checked by setproperty

if ~isnumeric(x1) | ~isnumeric(y1) | ~isnumeric(u) | ~isnumeric(v)
    error('Input of x1, y1, u and v should be numeric');
end

if min(size(x1)) == 1 & min(size(y1)) == 1 & min(size(u)) == 1 & min(size(v)) == 1
    % We are dealing with vectors
    if length(x1) == length(y1) & length(x1) == length(u) & length(x1) == length(v)
        % Vectors are consistent, make sure they have their length in the same direction:
        if size(x1,2)~=1
            x1=x1';
        end
        if size(y1,2)~=1
            y1=y1';
        end
        if size(u,2)~=1
            u=u';
        end
        if size(v,2)~=1
            v=v';
        end
    else
        error('Not all vectors have the same lengths, please check this');
    end
elseif min(size(x1)) > 1 & min(size(y1)) > 1 & min(size(u)) > 1 & min(size(v)) > 1
    if min(size(x1)) == min(size(y1)) & min(size(x1)) == min(size(u)) & min(size(x1)) == min(size(v))
        if max(size(x1)) == max(size(y1)) & max(size(x1)) == max(size(u)) & max(size(x1)) == max(size(v))
            %All matrices are of the same dimensions, now make sure the directions are as well..
            if find(max(size(x1))==size(x1))==2
                x1=x1';
            end
            if find(max(size(y1))==size(y1))==2
                y1=y1';
            end
            if find(max(size(u))==size(u))==2
                u=u';
            end
            if find(max(size(v))==size(v))==2
                v=v';
            end
            % Now lets turn them into vectors..
            x1 = x1(:);
            y1 = y1(:);
            u = u(:);
            v = v(:);
        else
            error('Input matrices are of different sizes');
        end
    else
        error('Input matrices are of different sizes');
    end
else
    error('Please check the input, it is not consistent (should be of the same size)');
end

% First convert x,y,u,v to x1,y1,x2,y2 using scale and/or uniformity
if strcmp(keywords.scaling_type,'uniform')
    x2 = x1+(keywords.uniform_length.*sin(atan2(u,v)));
    y2 = y1+(keywords.uniform_length.*cos(atan2(u,v)));
    if keywords.plot_to_current_figure == 1
        if isempty(cell2mat(strfind(varargin(cellfun(@ischar,varargin)),'plot_color')))
            varargin{1,end+1} = 'plot_color';
            varargin{1,end+1} = 'k';
            disp('Warning');
            disp(' ');
            disp('You are advised to specify ''plot_color'' yourself when using uniform scaling, it is now automatically set to black...');
        elseif strcmp(varargin{1,(find(strcmp(varargin(cellfun(@ischar,varargin)),'plot_color')==1)+1)},'magnitude')
            varargin{1,(find(cellfun(@isempty,strfind(varargin(cellfun(@ischar,varargin)),'plot_color'))==0)+1)} = 'k';
            disp('Warning');
            disp(' ');
            disp('When using uniform scaling you are not advised to use ''magnitude'' for the ''plot_color'' since the color of each');
            disp('arrow will depend on the relative error of the functions sin(e) and cos(ine), it is therefore automatically set to black, but you are free to change this...');
        end
    end
else % used for anything else then uniform (also magnitude, but doesn't matter
    x2 = x1+keywords.scaling.*u;
    y2 = y1+keywords.scaling.*v;
end

% Now call arrow_thick given the number of output arguments

if nargout < 2
    varargout = arrow_thick(x1,y1,x2,y2,varargin);
    % Please note that we get a struct back, not a (required) cell struct
    varargout = {varargout};
elseif nargout == 2
    [varargout{1} varargout{2}] = arrow_thick(x1,y1,x2,y2,varargin);
elseif nargout == 3
    [varargout{1} varargout{2} varargout{3}] = arrow_thick(x1,y1,x2,y2,varargin);
else
    error('More than 3 output arguments are not supported, please change this...');
end


end