% textsc() - places text in screen coordinates and places
%            a title at the top of the figure.
%
% Usage:
%   H = textsc(X,Y,TXT) places the text string, TXT
%   at the normalized coordinates X and Y.  H is the
%   handle to the text object.
%
%   H = textsc(TXT) or H = textsc(TXT,'title') or 
%   H = textsc('title', TXT) places a title at the top
%   of the figure window.  This is useful when you
%   want to place a single title above multiple
%   subplots.
%
% Notes: textsc creates an invisible AXES which occupies
% the entire FIGURE.  The units of the AXES are 
% normalized (range from 0 to 1).  textsc checks
% all the children of the current FIGURE to determine
% if an AXES exist which meets these criteria already
% exist.  If one does, then it places the text relative
% to that AXES.
%
% Author: John L. Galenski, January 21, 1994

% Written by John L. Galenski III
% All Rights Reserved  January 21, 1994
% LDM031695jlg

function H = textsc(x,y,txt)

% Basic error checking
if nargin < 2
    y = 'title';
end
if isstr(x) && strcmpi(x, 'title')
    x = y;
    y = 'title';
end

% Check to see if AXES already exist
ch = get(gcf,'Children');
if ~isempty(ch)
    try
        ind = cellfun(@(x)isequal('axes', x), get(ch, 'type'));
    catch
        ind = cellfun(@(x)isequal('axes', x), {get(ch, 'type')}); % fix Joe Dien bug 1538
    end
    if any(ind), ch = gca; end
end

ax = findobj(gcf,'Type','axes','Tag','TEXTSC');
if isempty(ax)
  ax = axes('Units','Normal','Position',[0 0 1 1], ...
            'Visible','Off','Tag','TEXTSC');
else
  axes(ax)
end

% Place the text
if ischar(y) && ischar(x) && strcmp(lower(y),'title')  % Subplot title
  txt = x;
  x = .5;
  tmp = text('Units','normal','String','tmp','Position',[0 0 0]);
  ext = get(tmp,'Extent');
  delete(tmp)
  H = ext(4);
  y = 1 - .60*H;
end
h = text(x,y,txt,'VerticalAlignment','Middle', ...
         'HorizontalAlignment','Center','interpreter', 'none' );

% Make the original AXES current
if ~isempty(ch)
  set(gcf,'CurrentAxes',ch(end));
end

% Check for output
if nargout == 1
  H = h;
end
