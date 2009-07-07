function [varargout] = plot_topo(dat, chanX, chanY, varargin)

% PLOT_TOPO calculates and plots topography (no layout)
%
% Use as
%   plot_topo(datavector, chanX, chanY, ...)
%   plot_topo(datavector, chanX, chanY, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'hpos'
%   'vpos'
%   'width'
%   'height'
%   'mask'
%   'shading'
%   'gridscale'

% Copyrights (C) 2009, Giovanni Piantoni
%
% $Log: not supported by cvs2svn $
% Revision 1.2  2009/06/02 15:36:25  giopia
% first implementation based on topoplot.m
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'gridscale', 'mask'});
hpos        = keyval('hpos', varargin);     if isempty(hpos);       hpos = 0;           end
vpos        = keyval('vpos', varargin);     if isempty(vpos);       vpos = 0;           end
width       = keyval('width', varargin);    if isempty(width);      width = 1;          end
height      = keyval('height', varargin);   if isempty(height);     height = 1;         end
gridscale   = keyval('gridscale', varargin);if isempty(gridscale);  gridscale = 67;     end; % 67 in original
shading     = keyval('shading', varargin);  if isempty(shading);    shading = 'flat';   end;
mask        = keyval('mask', varargin);

chanX = chanX * width  + hpos;
chanY = chanY * height + vpos;

hlim = [min(chanX) max(chanX)];
vlim = [min(chanY) max(chanY)];

xi         = linspace(hlim(1), hlim(2), gridscale);   % x-axis for interpolation (row vector)
yi         = linspace(vlim(1), vlim(2), gridscale);   % y-axis for interpolation (row vector)
[Xi,Yi,Zi] = griddata(chanX', chanY, dat, xi', yi, 'v4'); % Interpolate the topographic data

Zi(~mask) = NaN;

deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry
h = surface(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)), Zi, 'EdgeColor', 'none', 'FaceColor',shading);

% the (optional) output is the handle
if nargout == 1
  varargout{1} = h;
end
