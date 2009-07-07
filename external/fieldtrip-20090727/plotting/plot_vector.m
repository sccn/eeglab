function [varargout] = plot_vector(varargin)

% PLOT_VECTOR
%
% Use as
%   plot_vector(Y, ...)
%   plot_vector(X, Y, ...)
% where X and Y are similar as the input to the Matlab plot function.
%
% Additional options should be specified in key-value pairs and can be
%   'hpos'
%   'vpos'
%   'width'
%   'height'
%   'hlim'
%   'vlim'
%   'style'
%   'axis'    can be 'yes' or 'no'
%   'box'     can be 'yes' or 'no'
%
% Example use
%   plot_vector(randn(1,100), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0)

% Copyrights (C) 2009, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.8  2009/06/04 13:11:54  crimic
% added highlight option
%
% Revision 1.7  2009/06/02 15:42:52  giopia
% correct error in first if-statement and added varargout for handle
%
% Revision 1.6  2009/04/15 20:00:45  roboos
% small change in input parsing
%
% Revision 1.5  2009/04/14 19:48:28  roboos
% added keyvalcheck
%
% Revision 1.4  2009/04/14 18:55:51  roboos
% changed the handling of the input arguments for a closer resemblance to plot()
%
% Revision 1.3  2009/04/14 14:31:08  roboos
% many small changes to make it fully functional
%

if nargin>2 && all(cellfun(@isnumeric, varargin(1:2)))
  % the function was called like plot(x, y, ...)
  hdat = varargin{1};
  vdat = varargin{2};
  varargin = varargin(3:end);
else
  % the function was called like plot(y, ...)
  vdat = varargin{1};
  if any(size(vdat)==1)
    % ensure that it is a column vector
    vdat = vdat(:);
  end
  hdat = 1:size(vdat,1);
  varargin = varargin(2:end);
end

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'hlim', 'vlim', 'style', 'axis', 'box','highlight','highlightstyle'});
hpos   = keyval('hpos',   varargin);
vpos   = keyval('vpos',   varargin);
width  = keyval('width',  varargin);
height = keyval('height', varargin);
hlim   = keyval('hlim',   varargin); if isempty(hlim), hlim = 'maxmin'; end
vlim   = keyval('vlim',   varargin); if isempty(vlim), vlim = 'maxmin'; end
style  = keyval('style',  varargin); if isempty(style), style = '-'; end
axis   = keyval('axis',   varargin); if isempty(axis), axis = false; end
box    = keyval('box',    varargin); if isempty(box), box = false; end
highlight      = keyval('highlight',    varargin);
highlightstyle = keyval('highlightstyle',    varargin); if isempty(highlightstyle), highlightstyle = 'box'; end

% label  = keyval('label', varargin); % FIXME

if ischar(hlim)
  switch hlim
    case 'maxmin'
      hlim = [min(hdat) max(hdat)];
    case 'absmax'
      hlim = max(abs(hdat));
      hlim = [-hlim hlim];
    otherwise
      error('unsupported option for hlim')
  end % switch
end % if ischar

if ischar(vlim)
  switch vlim
    case 'maxmin'
      vlim = [min(vdat(:)) max(vdat(:))];
    case 'absmax'
      vlim = max(abs(vdat(:)));
      vlim = [-vlim vlim];
    otherwise
      error('unsupported option for vlim')
  end % switch
end % if ischar


if isempty(hpos) && ~isempty(hlim)
  hpos = (hlim(1)+hlim(2))/2;
end
if isempty(vpos) && ~isempty(vlim)
  vpos = (vlim(1)+vlim(2))/2;
end

if isempty(width) && ~isempty(hlim)
  width = hlim(2)-hlim(1);
end

if isempty(height) && ~isempty(vlim)
  height = vlim(2)-vlim(1);
end


% first shift the horizontal axis to zero
hdat = hdat - (hlim(1)+hlim(2))/2;
% then scale to length 1
hdat = hdat ./ (hlim(2)-hlim(1));
% then scale to the new width
hdat = hdat .* width;
% then shift to the new horizontal position
hdat = hdat + hpos;
% first shift the vertical axis to zero
vdat = vdat - (vlim(1)+vlim(2))/2;
% then scale to length 1
vdat = vdat ./ (vlim(2)-vlim(1));
% then scale to the new width
vdat = vdat .* height;
% then shift to the new vertical position
vdat = vdat + vpos;

if ~isempty(highlight)
  switch highlightstyle
    case 'box'
      % find the sample number where the highligh begins and ends
      if ~islogical(highlight)
        highlight=logical(highlight);
        warning('converting mask to logical values')
      end
      begsample = find(diff([0 highlight 0])== 1);
      endsample = find(diff([0 highlight 0])==-1)-1;
      for i=1:length(begsample)
        begx = hdat(begsample(i));
        endx = hdat(endsample(i));
        plot_box([begx endx vpos-height/2 vpos+height/2], 'facecolor', [.6 .6 .6], 'edgecolor', 'none');
      end
    case 'thickness'
      error('unsupported highlightstyle')
    case 'opacity'
      error('unsupported highlightstyle')
    otherwise
      error('unsupported highlightstyle')
  end % switch highlightstyle
end

h = plot(hdat, vdat, style);



if istrue(box)
  boxposition = zeros(1,4);
  % this plots a box around the original hpos/vpos with appropriate width/height
  boxposition(1) = hpos - width/2;
  boxposition(2) = hpos + width/2;
  boxposition(3) = vpos - height/2;
  boxposition(4) = vpos + height/2;
  plot_box(boxposition, 'facecolor', 'none', 'edgecolor', 'k');
  
  % this plots a box around the complete data
  % boxposition(1) = hlim(1);
  % boxposition(2) = hlim(2);
  % boxposition(3) = vlim(1);
  % boxposition(4) = vlim(2);
  % plot_box(boxposition, 'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim);
end

if istrue(axis)
  %   X = hlim;
  %   Y = [0 0];
  %   plot_line(X, Y, 'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim);
  %   str = sprintf('%g', hlim(1)); plot_text(X(1), Y(1), str, 'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim);
  %   str = sprintf('%g', hlim(2)); plot_text(X(2), Y(2), str, 'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim);
  
  X = [hpos-width/2  hpos+width/2];
  Y = [vpos vpos];
  plot_line(X, Y);
  str = sprintf('%g', hlim(1)); plot_text(X(1), Y(1), str);
  str = sprintf('%g', hlim(2)); plot_text(X(2), Y(2), str);
  
  X = [hpos hpos];
  Y = [vpos-height/2 vpos+height/2];
  plot_line(X, Y);
  str = sprintf('%g', vlim(1)); plot_text(X(1), Y(1), str);
  str = sprintf('%g', vlim(2)); plot_text(X(2), Y(2), str);
end



% the (optional) output is the handle
if nargout == 1;
  varargout{1} = h;
end

