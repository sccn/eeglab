function [varargout] = plot_text(X, Y, str, varargin)

% PLOT_TEXT

% Copyrights (C) 2009, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.4  2009/06/02 15:40:36  giopia
% added varargout to pass handle
%
% Revision 1.3  2009/04/14 19:48:28  roboos
% added keyvalcheck
%
% Revision 1.2  2009/04/14 14:31:08  roboos
% many small changes to make it fully functional
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'hlim', 'vlim', 'color', 'fontsize', 'fontname'});
hpos        = keyval('hpos',      varargin);
vpos        = keyval('vpos',      varargin);
width       = keyval('width',     varargin);
height      = keyval('height',    varargin);
hlim        = keyval('hlim',      varargin);
vlim        = keyval('vlim',      varargin);
color       = keyval('color',     varargin);  if isempty(color), color = 'k'; end
fontsize    = keyval('fontsize',  varargin);
fontname    = keyval('fontname',  varargin);

abc = axis;
if isempty(hlim)
  hlim = abc([1 2]);
end

if isempty(vlim)
  vlim = abc([3 4]);
end

if isempty(hpos);
  hpos = (hlim(1)+hlim(2))/2;
end

if isempty(vpos);
  vpos = (vlim(1)+vlim(2))/2;
end

if isempty(width),
  width = hlim(2)-hlim(1);
end

if isempty(height),
  height = vlim(2)-vlim(1);
end

% first shift the horizontal axis to zero
X = X - (hlim(1)+hlim(2))/2;
% then scale to length 1
X = X ./ (hlim(2)-hlim(1));
% then scale to the new width
X = X .* width;
% then shift to the new horizontal position
X = X + hpos;

% first shift the vertical axis to zero
Y = Y - (vlim(1)+vlim(2))/2;
% then scale to length 1
Y = Y ./ (vlim(2)-vlim(1));
% then scale to the new width
Y = Y .* height;
% then shift to the new vertical position
Y = Y + vpos;

h = text(X, Y, str);
set(h, 'HorizontalAlignment', 'center');
% set(h, 'VerticalAlignment', 'middle'); % this is already the default
set(h, 'Color', color);
if ~isempty(fontsize), set(h, 'FontSize', fontsize); end
if ~isempty(fontname), set(h, 'FontName', fontname); end

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = h;
end
