function varargout = arrow_thick(x1,y1,x2,y2,varargin)
% arrow_thick plot fancy arrows like arrow and quiver
%__________________________________________________________________________
%Introduction
%
%ARROW_THICK allows users to plot arrows just like the syntax of the
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
%arrow_thick uses start (x1,y1) and end position (x2,y2) to generate 2D
%plots of arrows as filled patches. Thus, the function works just like the
%'arrow' function. The start & end (input) locations can be specified on 
%an (irregular) grid/mesh or basic matlab vectors. When using (X,Y) along
%with (U,V) vector components, one is referred to quiver_thick, which works
%just like the well known function 'quiver' with (X,Y,U,V). By default,
%arrows are colored (filled) according to their magnitude.
%
%Syntax:
%
%<output> = arrow_thick(x1,y1,x2,y2,<keyword,value>);
%
%Apart from x1, y1, x2 & y2 please note that all other in- and output
%variables (depicted by <...>) are optional
%
%Input variables:
%
%x1                        Required: Start X-coordinates of arrows, can be
%                          any vector or matrix [MxN] M,N = 1:any integer
%y1                        Required: Start Y-coordinates of arrows, can be
%                          any vector or matrix [MxN] M,N = 1:any integer
%x2                        Required: End X-coordinates of arrows, can be
%                          any vector or matrix [MxN] M,N = 1:any integer
%y2                        Required: End Y-coordinates of arrows, can be
%                          any vector or matrix [MxN] M,N = 1:any integer
%
%All <keyword,value> pairs below are optional and shown as <name,default>
%
%Simply call arrow_thick to get a list of optional keywords and their defaults
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
%original start and end locations
%__________________________________________________________________________
%Examples
%_________
%Example 1 - All default settings:
%
%arrow_thick(zeros(1,36),...
%            zeros(1,36),...
%            (cosd(1:10:360)).*(36:-1:1),...
%            (sind(1:10:360)).*(36:-1:1));
%box on; axis tight; set(gca,'XTick',[],'YTick',[]);
%_________
%Example 2 - Making normal arrow/quiver plots
%
%x=0:10; y=0:10; [X1,Y1]=meshgrid(x,y); X2=X1+rand(11,11); Y2=Y1+2*rand(11,11);
%[X_data,Y_data] = arrow_thick(X1,Y1,X2,Y2,...
%                 'arrow_thickness',0,...
%                 'arrowhead_length_option','absolute',...
%                 'plot_color','k');
%box on; axis tight; set(gca,'XTick',[],'YTick',[]);
%_________
%Example 3 - Getting the arrow data and plot afterwards
%
%figure; plot(0:10,0:10,'k--'); hold on; plot(0:10,10:-1:0,'k--');
%OutputData = arrow_thick([6 6 4 4],[6 4 4 6],[9 9 1 1],[9 1 1 9],...
%                         'plot_to_current_figure',0);
%p1 = patch(OutputData.arrows_x_values,OutputData.arrows_y_values,'r');
%box on; axis equal; axis tight;
%set(p1,'FaceAlpha',0.6); set(gca,'XTick',[],'YTick',[]);
%_________
%Example 4 - Using many options
%
%X1=rand(999,1)*20-10; Y1=rand(999,1)*20-10; X2=X1+0.2*X1; Y2=Y1+0.2*Y1;
%text(0,0,'KABOOM!','horizontalalignment','center',...
%         'fontsize',40,'fontweight','bold','color','r'); hold on;
%OutputData = arrow_thick(X1,Y1,X2,Y2,...
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
%See also: quiver_thick patch arrow quiver curvec

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

% Default keywords from arrow_thick:

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

% Default keywords from quiver_thick
keywords.scaling                 = 1;
keywords.scaling_type            = 'magnitude';
keywords.uniform_length          = 1;

if nargin == 0
    default_keywords = keywords
    disp('Please note that the keywords ''scaling'', ''scaling_type'' & ''uniform_length'' are only used when calling quiver_thick');
    return
end

% If a call comes in from quiver_thick, varargin is a 1 by 1 cell:
if size(varargin,1)==1 & size(varargin,2)==1
    varargin = varargin{:};
    if odd(length(varargin))
        error('Odd number of input parameters for <keyword,value> pairs supplied, please check this');
    end
elseif odd(length(varargin)) % else we just check for consistency
    error('Odd number of input parameters for <keyword,value> pairs supplied, please check this');
end

% Set user defined changes to default keywords:

keywords = setproperty(keywords,varargin);

% Do some basic checks for input, all <keyword,value> pairs are checked by setproperty

if ~isnumeric(x1) | ~isnumeric(y1) | ~isnumeric(x2) | ~isnumeric(y2)
    error('Input of x1, y1, x2 and y2 should be numeric');
end

if min(size(x1)) == 1 & min(size(y1)) == 1 & min(size(x2)) == 1 & min(size(y2)) == 1
    % We are dealing with vectors
    if length(x1) == length(y1) & length(x1) == length(x2) & length(x1) == length(y2)
        % Vectors are consistent, make sure they have their length in the same direction:
        if size(x1,2)~=1
            x1=x1';
        end
        if size(y1,2)~=1
            y1=y1';
        end
        if size(x2,2)~=1
            x2=x2';
        end
        if size(y2,2)~=1
            y2=y2';
        end
    else
        error('Not all vectors have the same lengths, please check this');
    end
elseif min(size(x1)) > 1 & min(size(y1)) > 1 & min(size(x2)) > 1 & min(size(y2)) > 1
    if min(size(x1)) == min(size(y1)) & min(size(x1)) == min(size(x2)) & min(size(x1)) == min(size(y2))
        if max(size(x1)) == max(size(y1)) & max(size(x1)) == max(size(x2)) & max(size(x1)) == max(size(y2))
            %All matrices are of the same dimensions, now make sure the directions are as well..
            if find(max(size(x1))==size(x1))==2
                x1=x1';
            end
            if find(max(size(y1))==size(y1))==2
                y1=y1';
            end
            if find(max(size(x2))==size(x2))==2
                x2=x2';
            end
            if find(max(size(y2))==size(y2))==2
                y2=y2';
            end
            % Now lets turn them into vectors..
            x1 = x1(:);
            y1 = y1(:);
            x2 = x2(:);
            y2 = y2(:);
        else
            error('Input matrices are of different sizes');
        end
    else
        error('Input matrices are of different sizes');
    end
else
    error('Please check the input, it is not consistent (should be of the same size)');
end

% get some basic arrow stats

arrow_lengths = sqrt((diff([x1 x2],1,2).^2)+(diff([y1 y2],1,2).^2));
arrow_angle   = atan2(diff([x1 x2],1,2),diff([y1 y2],1,2));
arrow_x_cos   = cos(arrow_angle);
arrow_y_sin   = sin(arrow_angle);

% Convert to absolute arrow width:

if strcmp(keywords.arrow_thickness_option,'relative')
    abs_arrow_widths = keywords.arrow_thickness.*arrow_lengths;
elseif strcmp(keywords.arrow_thickness_option,'absolute')
    abs_arrow_widths = repmat(keywords.arrow_thickness,size(x1,1),size(x1,2));
else
    disp(['Option ' keywords.arrow_thickness_option ' for keyword ''arrow_thickness_option'' is unknown, using default (''relative'')']);
    abs_arrow_widths = keywords.arrow_thickness.*arrow_lengths;
end

% Convert to absolute arrowhead length:

if strcmp(keywords.arrowhead_length_option,'relative')
    abs_arrowhead_length = keywords.arrowhead_length.*arrow_lengths;
elseif strcmp(keywords.arrowhead_length_option,'absolute')
    abs_arrowhead_length = repmat(keywords.arrowhead_length,size(x1,1),size(x1,2));
else
    disp(['Option ' keywords.arrowhead_length_option ' for keyword ''arrowhead_length_option'' is unknown, using default (''relative'')']);
    abs_arrowhead_length = keywords.arrowhead_length.*arrow_lengths;
end

% Do some checks (details behind call)

abs_arrow_widths     = max([zeros(size(x1,1),size(x1,2)) abs_arrow_widths],[],2); % makes sure the width is at least 0
abs_arrowhead_length = max([zeros(size(x1,1),size(x1,2)) min([0.9*arrow_lengths abs_arrowhead_length],[],2)],[],2); % makes sure the arrowhead length is at max 0.9 times the length of the arrow itself and at least zero
abs_head_widths      = max([abs_arrow_widths abs_arrowhead_length.*tand(keywords.max_arrowhead_angle/2)],[],2); % makes sure the arrowhead is at least at wide as the width of the arrow

% Generate all arrows in 2 calls (used for patch plotting)

output.arrows_x_values = [x1+(abs_arrow_widths.*arrow_x_cos) x1 x1-(abs_arrow_widths.*arrow_x_cos) x2-(arrow_y_sin.*abs_arrowhead_length)-(abs_arrow_widths.*arrow_x_cos) x2-(arrow_y_sin.*abs_arrowhead_length)-(abs_head_widths.*arrow_x_cos) x2 x2-(arrow_y_sin.*abs_arrowhead_length)+(abs_head_widths.*arrow_x_cos) x2-(arrow_y_sin.*abs_arrowhead_length)+(abs_arrow_widths.*arrow_x_cos) x1+(abs_arrow_widths.*arrow_x_cos)]';
output.arrows_y_values = [y1-(abs_arrow_widths.*arrow_y_sin) y1 y1+(abs_arrow_widths.*arrow_y_sin) y2-(arrow_x_cos.*abs_arrowhead_length)+(abs_arrow_widths.*arrow_y_sin) y2-(arrow_x_cos.*abs_arrowhead_length)+(abs_head_widths.*arrow_y_sin) y2 y2-(arrow_x_cos.*abs_arrowhead_length)-(abs_head_widths.*arrow_y_sin) y2-(arrow_x_cos.*abs_arrowhead_length)-(abs_arrow_widths.*arrow_y_sin) y1-(abs_arrow_widths.*arrow_y_sin)]';

% Plot to the current figure using hold on and axis equal (if not turned off by the user)

if keywords.plot_to_current_figure
    hold on;
    if strcmp(keywords.plot_color,'magnitude')
        output.arrows_handle = patch(output.arrows_x_values,output.arrows_y_values,arrow_lengths');
        colormap(keywords.colormap)
        if keywords.plot_colorbar == 1
            output.colorbar_handle = colorbar;
        end
    else
        output.arrows_handle = patch(output.arrows_x_values,output.arrows_y_values,keywords.plot_color);
    end
    if keywords.axis_equal
        axis equal
    end
    
    set(output.arrows_handle,'EdgeColor',keywords.EdgeColor);
    set(output.arrows_handle,'EdgeAlpha',keywords.EdgeAlpha);
    if keywords.FaceAlpha~=1
        set(output.arrows_handle,'FaceAlpha',keywords.FaceAlpha);
    end
    if keywords.LineWidth~=0.5
        set(output.arrows_handle,'LineWidth',keywords.LineWidth);
    end
    if strcmp(keywords.LineStyle,'-')==0
        set(output.arrows_handle,'LineStyle',keywords.LineStyle);
    end
end

% Generate output variables based on number of output variables specified by the user

if nargout == 0
    varargout{1} = output;
elseif nargout == 1
    varargout{1} = output;
elseif nargout == 2
    varargout{1} = output.arrows_x_values;
    varargout{2} = output.arrows_y_values;
else
    for ii=1:nargout
        if ii==1
            varargout{ii} = output.arrows_x_values;
        elseif ii==2
            varargout{ii} = output.arrows_y_values;
        elseif ii==3
            if keywords.plot_to_current_figure
                varargout{ii} = output.arrows_handle;
            else
                varargout{ii} = ['No handle assigned since plotting is turned off'];
            end
        elseif ii==4
            if keywords.plot_colorbar == 1
                if keywords.plot_to_current_figure
                    varargout{ii} = output.colorbar_handle;
                else
                    varargout{ii} = ['No handle assigned since plotting is turned off'];
                end
            else
                varargout{ii} = ['No colorbar handle was made as ''plot_colorbar'' = 0'];
            end
        else
            varargout{ii} = ['Output argument no. ' num2str(ii) ' is not assigned, too many'];
        end
    end
end

end