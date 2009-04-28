function [cmap] = ColorMapsShow

% colorMapsShow - display/select available Matlab colormaps
%
% Example:  figure, colormap(ColorMapsShow), colorbar
%
% See also colorMapsMake, colorbar
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:50 $

% Licence:  GNU GPL, no implied or express warranties
% History:  01/99 Uilke Stelwagen, Copyright (C) 1992-1999
%           Institute of Applied Physics, TNO-TPD, The Netherlands.
%           09/01 Darren.Weber_at_radiology.ucsf.edu
%                 - downloaded function from the mathworks, distrib.
%                   under GPL
%                 - adapted function to gui selection of colormap
%                   and return of a prefered map in cmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maps = [
   'hsv        - hue-saturation-value           '
   'hot        - black-red-yellow-white         '
   'gray       - linear gray-scale              '
   'bone       - gray-scale with blue tinge     '
   'copper     - linear copper-tone             '
   'pink       - pastel pinks                   '
   'white      - all white                      '
   'flag       - red, white, blue, & black      '
   'lines      - color map of the line colors   '
   'colorcube  - enhanced color-cube            '
   'jet        - variant HSV (Matlab default)   '
   'prism      - prism                          '
   'cool       - cyan and magenta               '
   'autumn     - red and yellow                 '
   'spring     - magenta and yellow             '
   'winter     - blue and green                 '
   'summer     - green and yellow               '
];

[map,ok] = listdlg('PromptString','Double-Click ColorMap',...
                   'SelectionMode','single',...
                   'ListString',maps,...
                   'CancelString','Done',...
                   'uh',30,'fus',5,'ffs',5,...
                   'ListSize',[200 330]);

figure('menubar','none','NumberTitle','off',...
       'NextPlot','replace');

h = subplot(143);
colorbar(h);

while ok
    
    colormap(deblank(maps(map,1:10)));
    %title(   deblank(maps(map,:)));
    
    [map,ok] = listdlg('PromptString','Double-Click ColorMap',...
                       'SelectionMode','single',...
                       'InitialValue',map,...
                       'ListString',maps,...
                       'CancelString','Done',...
                       'uh',30,'fus',5,'ffs',5,...
                       'ListSize',[200 330]);
                   
end
cmap = colormap;
close gcf
