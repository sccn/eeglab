function [varargout] = componentbrowser(cfg, comp)

% COMPONENTBROWSER plots topography and activations of ICA components
%
% Use as
%   componentbrowser(cfg, comp)
% where comp is a FieldTrip structure obtained from COMPONENTANALYSIS.
%
% The configuration has the following parameters:
% cfg.layout = layout from PREPARE_LAYOUT (required)
% cfg.comp   = a vector with the components to plot (ex. 1:10) (optional)
% cfg.trial  = choose which trial to plot first (optional, only one trial)

% Copyright (C) 2009, Giovanni Piantoni
%
% $Log: not supported by cvs2svn $
% Revision 1.3  2009/06/19 15:11:00  giopia
% allows scroll through components
%
% Revision 1.2  2009/06/03 14:00:26  roboos
% fixed cfg.lay, should be cfg.layout
%
% Revision 1.1  2009/06/02 15:48:58  giopia
% first implementation, plot topoplot, activations and simple interactive
%

fieldtripdefs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that the data comes from componentanalysis
comp = checkdata(comp, 'datatype', 'comp');

% set the defaults:
if ~isfield(cfg, 'comp'),   cfg.comp  = 1:10; end
if ~isfield(cfg, 'trial'),  cfg.trial = 1;    end

if numel(cfg.trial) > 1,
  warning('componentbrowser:cfg_onetrial', 'only one trial can be plotted at the time');
  cfg.trial = cfg.trial(1);
end

% Read or create the layout that will be used for plotting:
[cfg.layout] = prepare_layout(cfg, comp);

% Identify the channels to plot
[labels, cfg.chanidx] = intersect(comp.topolabel, cfg.layout.label); % in case channels are missing
if isempty(cfg.chanidx)
  error('componentbrowser:labelmismatch', 'The channel labels in the data do not match the labels of the layout');
end

% fixed variables
cfg.shift   = 1.2;   % distance between topoplots
cfg.gridscale   = 67;    % default parameter from topoplot
[cfg.mask] = createmask(cfg.layout, cfg.gridscale);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure and assign userdata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create figure
[cfg.h] = figure('uni','pix', 'name', 'componentbrowser', 'vis', 'off', 'numbertitle', 'off');
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buttons and Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scroll components
uicontrol(cfg.h,'uni','pix','pos',[105 5 25 18],'str','-',...
  'call',{@plottopography, comp});

cfg.ncomp = uicontrol(cfg.h,'sty','text','uni','pix','pos',[130 5 150 18],...
  'str',['comp n.' num2str(cfg.comp(1)) '-' num2str(cfg.comp(end))]);

uicontrol(cfg.h,'uni','pix','pos',[280 5 25 18],'str','+',...
  'call',{@plottopography, comp});

% scroll trials
uicontrol(cfg.h,'uni','pix','pos',[355 5 25 18],'str','<<',...
  'call',{@plotactivation, comp});

cfg.ntrl = uicontrol(cfg.h,'sty','text','uni','pix','pos',[380 5 70 18],...
  'str',['trial n.' num2str(cfg.trial)]);

uicontrol(cfg.h,'uni','pix','pos',[450 5 25 18],'str','>>',...
  'call',{@plotactivation, comp});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First callback and final adjustments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first call of the two plotting functions
plottopography([], cfg, comp)
plotactivation([], cfg, comp)

% final adjustments
set(cfg.h, 'vis', 'on')
axis equal
axis off
hold off

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = cfg.h;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTOPOGRAPHY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plottopography(h, cfg, comp)
% now plottopography is not associated with a callback, but it might in
% the future

if isempty(h) % when called in isolation
  set(cfg.h, 'user', cfg)
else
  cfg = get(get(h, 'par'), 'user');

  % which button has been pressed
  if intersect(h, findobj(cfg.h, 'str', '+'))

    cfg.comp = cfg.comp + numel(cfg.comp);
    if cfg.comp(end) > size(comp.label,1)
      cfg.comp = cfg.comp - (cfg.comp(end) - size(comp.label,1));
    end

  elseif intersect(h, findobj(cfg.h, 'str', '-'))

    cfg.comp = cfg.comp - numel(cfg.comp);
    if cfg.comp(1) < 1
      cfg.comp = cfg.comp - cfg.comp(1) + 1;
    end

  end
end

set(cfg.ncomp, 'str', ['comp n.' num2str(cfg.comp(1)) '-' num2str(cfg.comp(end))])
drawnow
delete(findobj(cfg.h, 'tag', 'comptopo'))

cnt = 0;
for k = cfg.comp
  cnt = cnt + 1;

  % write number of the component on the left
  h_text(cnt) = plot_text(-2.5, -cnt*cfg.shift, ['n. ' num2str(cfg.comp(cnt))]);

  % plot only topography (no layout)
 h_topo(cnt) = plot_topo(comp.topo(cfg.chanidx, k), cfg.layout.pos(cfg.chanidx,1), cfg.layout.pos(cfg.chanidx,2), ...
    'hpos', -1, 'vpos', -cnt*cfg.shift, 'mask', cfg.mask, 'gridscale', cfg.gridscale);
  % plot layout
  plot_lay(cfg.layout, 'hpos', -1, 'vpos', -cnt*cfg.shift, 'point', false, 'box', false, 'label', false, 'mask', true, 'verbose', false);
end

set(h_text, 'tag', 'comptopo')
set(h_topo, 'tag', 'comptopo')

plotactivation([], cfg, comp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTACTIVATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotactivation(h, cfg, comp)
% plotactivation can be called in isolation or by buttondownfcn
% cfg is stored in 'user' of the main figure

if isempty(h) % when called in isolation
  set(cfg.h, 'user', cfg)
else
  cfg = get(get(h, 'par'), 'user');

  % which button has been pressed
  if intersect(h, findobj(cfg.h, 'str', '>>'))

    cfg.trial = cfg.trial + 1;
    if cfg.trial > size(comp.trial,2)
      cfg.trial = size(comp.trial,2);
    end

  elseif intersect(h, findobj(cfg.h, 'str', '<<'))

    cfg.trial = cfg.trial - 1;
    if cfg.trial < 1
      cfg.trial = 1;
    end

  end
end

set(cfg.ntrl,'str',['trial n. ' num2str(cfg.trial)])
drawnow
delete(findobj(cfg.h,'tag', 'activations'));

hold on
cnt = 0;
for k = cfg.comp
  cnt = cnt + 1;

  % plot the activations
  h_act(cnt) = plot_vector(comp.trial{cfg.trial}(k,:), 'hpos', 6 , 'vpos', -cnt*cfg.shift, 'width', 12, 'height', 1, 'box', true);
end

h_inv = plot(6+12+1, -cnt*cfg.shift, '.'); %
set(h_inv, 'vis', 'off')

set(h_act, 'tag', 'activations')
set(cfg.h, 'user', cfg)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATEMASK: create anatomical mask, only one for all the topoplots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask] = createmask(lay, gridscale)
% calculate anatomical mask only once, based only on lay

% find limits for interpolation:
xmin = +inf;
xmax = -inf;
ymin = +inf;
ymax = -inf;
for i=1:length(lay.mask)
  xmin = min([xmin; lay.mask{i}(:,1)]);
  xmax = max([xmax; lay.mask{i}(:,1)]);
  ymin = min([ymin; lay.mask{i}(:,2)]);
  ymax = max([ymax; lay.mask{i}(:,2)]);
end

xi = linspace(xmin, xmax, gridscale);   % x-axis for interpolation (row vector)
yi = linspace(ymin, ymax, gridscale);   % y-axis for interpolation (row vector)
Xi =  ones(gridscale,1)*xi;
Yi = (ones(gridscale,1)*yi)';

% apply anatomical mask to the data, i.e. that determines that the interpolated data outside the circle is not displayed
mask = false(gridscale);
for i=1:length(lay.mask)
  lay.mask{i}(end+1,:) = lay.mask{i}(1,:); % force them to be closed
  mask(inside_contour([Xi(:) Yi(:)], lay.mask{i})) = true;
end


