function label = select_channel(h, eventdata, varargin)

% SELECT_CHANNEL is a helper function that can be used as callback function
% in a figure. It allows the user to select a channel. The channel labels
% are returned.
%
% Use as
%   label = select_channel(h, eventdata, ...)
% The first two arguments are automatically passed by Matlab to any
% callback function. Additional options should be specified in key-value
% pairs and can be
%   'callback'  = function handle to be executed after channels have been
%                 selected
% You can pass additional arguments to the callback function in a cell-array
% like {@function_handle,arg1,arg2}
%
% Example
%   % create a figure
%   lay = prepare_layout([])
%   plot_lay(lay)
%
%   % add the required guidata
%   info       = guidata(gcf)
%   info.x     = lay.pos(:,1);
%   info.y     = lay.pos(:,2);
%   info.label = lay.label
%   guidata(gcf, info)
%
%   % add this function ass callback
%   set(gcf, 'WindowButtonUpFcn', {@select_channel, 'callback', @disp})
%
% Subsequently you can click in the figure and you'll see that the disp
% function is executed as callback and that it displays the selected
% channels.

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.2  2009/05/12 18:10:43  roboos
% added handling of follow-up callback function
%
% Revision 1.1  2009/05/12 12:49:33  roboos
% new implementation of helper function that can be used as callback in a figure
%

% get optional input arguments
callback   = keyval('callback', varargin{:});

pos = get(gca, 'CurrentPoint');
x   = pos(1,1);
y   = pos(1,2);

info     = guidata(h);
chan_x   = info.x;
chan_y   = info.y;
chan_lab = info.label;

% compute the distance between the clicked point and all channels
chan_pos  = [chan_x chan_y];
chan_dist = dist(chan_pos');
chan_dist = triu(chan_dist, 1);
chan_dist = chan_dist(:);
chan_dist = chan_dist(chan_dist>0);
% allow for some tolerance in the clicking
chan_dist = median(chan_dist);
tolerance = 0.3*chan_dist;

dx = chan_x - x;
dy = chan_y - y;
dd = sqrt(dx.^2 + dy.^2);
[d, i] = min(dd);
if d<tolerance
  label = chan_lab{i};
  fprintf('channel "%s" selected\n', label);
else
  label = {};
  fprintf('no channel selected\n');
end

if ~isempty(callback)
  if isa(callback, 'cell')
    % the callback specifies a function and additional arguments
    funhandle = callback{1};
    funargs   = callback(2:end);
    feval(funhandle, label, funargs{:});
  else
    % the callback only specifies a function
    funhandle = callback;
    feval(funhandle, label);
  end
end

if 0
  cfg = info.cfg;
  cfg.cohrefchannel = label;
  data = info.data;
  figure; multiplotER(cfg, data);
end

uiresume; % ??
