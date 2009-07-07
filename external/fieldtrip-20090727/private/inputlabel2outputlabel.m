function [outputlabel, outputindex] = inputlabel2outputlabel(cfg, freq)

% INPUTLABEL2OUTPUTLABEL is a subfunction which outputs the cell-arrays
% outputlabel and the corresponding outputindex, and defines how the 
% channels in the original data have to be combined, to provide the 
% wished for combination of the channels, as defined in cfg.combinechan
%
% Configuration-options are:
%   cfg.combinechan = 'planar' combines the horizontal and vertical planar-gradients
%                     'pseudomeg' one gradiometer versus the rest
%   TODO: more flexible way of combining, e.g. by providing a cell-array 

% $Log: not supported by cvs2svn $
% Revision 1.2  2006/06/23 10:51:02  jansch
% changed format of outputlabel
%
% Revision 1.1  2005/08/15 15:16:19  jansch
% First implementation. Moved out of freqdescriptives.m
%

if ~isfield(cfg, 'combinechan'), cfg.combinechan = 'no'; end;

if strcmp(cfg.combinechan, 'no')
  % the output labels are similar to the input labels
  outputlabel = {};
  outputindex = {};
  for i=1:length(freq.label)
    outputindex{i} = i;
    outputlabel(i) = freq.label(i);
  end
elseif strcmp(cfg.combinechan, 'planar')
  % find the combination of horizontal and vertical channels that should be combined
  planar      = planarchannelset(freq);
  [sel_dH, H] = match_str(freq.label, planar(:,1));  % indices of the horizontal channels
  [sel_dV, V] = match_str(freq.label, planar(:,2));  % indices of the vertical   channels
  lab_dH      = freq.label(sel_dH);
  lab_dV      = freq.label(sel_dV);

  if length(sel_dH)~=length(sel_dV)
    error('not all planar channel combinations are complete')
  end
  
  % find the other channels that are present in the data
  sel_other = setdiff(1:length(freq.label), [sel_dH(:)' sel_dV(:)']);
  lab_other = freq.label(sel_other);
  
  % define the channel names after combining the planar combinations
  % they should be sorted according to the order of the horizontal planar channels in the data
  [dum, sel_planar] = match_str(freq.label, planar(:,1));
  lab_comb          = planar(sel_planar,3);
  [dum, sel_comb]   = match_str(planar(H,3),planar(V,3));
  
  outputlabel = {};
  outputindex = {};

  outputlabel = [outputlabel; lab_other];
  outputindex = [outputindex; mat2cell(sel_other', ones(length(sel_other),1), 1)];

  outputlabel = [outputlabel; lab_comb];
  outputindex = [outputindex; mat2cell([sel_dH sel_dV(sel_comb)], ones(length(sel_dH),1), 2)];
  
elseif strcmp(cfg.combinechan, 'pseudomeg'),
  outputlabel = {};
  outputindex = {};
  for i=1:length(freq.label)
    outputlabel{i} = ['all-',freq.label{i}];
    outputindex{i} = [setdiff(1:length(freq.label),i)];
  end
  for i=1:length(freq.label)
    outputlabel{end+1} = freq.label{i};
    outputindex{end+1} = i;
  end
elseif iscell(cfg.combinechan(1)),
  for i=1:length(cfg.combinechan)
    outputindex{i} = match_str(freq.label, cfg.combinechan{i});
    outputindex{i} = outputindex{i}(:)';
    outputlabel{i} = cell2mat(freq.label(outputindex{i})');
  end
else
  error('unknown combination method');
end
