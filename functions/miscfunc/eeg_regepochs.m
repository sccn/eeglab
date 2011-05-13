% eeg_regepochs() - Convert a continuous dataset into consecutive epochs of 
%                   a specified regular length by adding dummy events of type 
%                   and epoch the data around these events. Alternatively
%                   only insert events for extracting these epochs.
%                   May be used to split up continuous data for
%                   artifact rejection followed by ICA decomposition.
%                   The computed EEG.icaweights and EEG.icasphere matrices
%                   may then be exported to the continuous dataset and/or 
%                   to its epoched descendents.
% Usage:
%     >> EEGout = eeg_regepochs(EEG); % use defaults
%     >> EEGout = eeg_regepochs(EEG, 'key', val, ...); 
%
% Required input:
%     EEG                 - EEG continuous data structure (EEG.trials = 1)
%
% Optional inputs:
%     'recurrence' - [in s] the regular recurrence interval of the dummy
%                    events used as time-locking events for the 
%                    consecutive epochs {default: 1 s}
%     'limits'     - [minsec maxsec] latencies relative to the time-locking
%                    events to use as epoch boundaries. Stated epoch length 
%                    will be reduced by one data point to avoid overlaps 
%                    {default: [0 recurrence_interval]}
%     'rmbase'     - [NaN|latency] remove baseline. NaN does not remove
%                    baseline. 0 remove the pre-0 baseline. To
%                    remove the full epoch baseline, enter a value
%                    larger than the upper epoch limit. Default is 0.
%     'eventtype'  - [string] name for the event type. Default is 'X'
%     'extractepochs' - ['on'|'off'] extract data epochs ('on') or simply
%                    insert events ('off'). Default is 'on'.
%
% Outputs:
%     EEGout              - the input EEG structure separated into consecutive 
%                           epochs.
%
% See also: pop_editeventvals(), pop_epoch(), rmbase();
%
% Authors: Hilit Serby, Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, Sep 02, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, Sep 02, 2005, hilit@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function EEG = eeg_regepochs(EEG, varargin)

if nargin < 1
    help eeg_regepochs;
    return;
end;

if length(EEG) > 1
    EEG = pop_mergeset(EEG, [1:length(EEG)]);
end;

% test input variables
% --------------------
if ~isstruct(EEG) | ~isfield(EEG,'event')
   error('first argument must be an EEG structure')
elseif EEG.trials > 1
   error('input dataset must be continuous data (1 epoch)');
end

if nargin > 1 && ~isstr(varargin{1})
    options = {};
    if nargin >= 2, options = { options{:} 'recurrence' varargin{1} }; end;
    if nargin >= 3, options = { options{:} 'limits'     varargin{2} }; end;
    if nargin >= 4, options = { options{:} 'rmbase'     varargin{3} }; end;
else
    options = varargin;
end;
g = finputcheck(options, { 'recurrence'    'real'  []  1;
                            'limits'        'real'  []  [0 1];
                            'rmbase'        'real'  []  0;
                            'eventtype'     'string' {} 'X';
                            'extractepochs' 'string' { 'on','off' } 'on' }, 'eeg_regepochs');
if isstr(g), error(g); end;

if g.recurrence < 0 | g.recurrence > EEG.xmax
  error('recurrence_interval out of bounds');
end

if nargin < 3
  g.limits = [0 g.recurrence];
end

if length(g.limits) ~= 2 | g.limits(2) <= g.limits(1) 
   error('epoch limits must be a 2-vector [minsec maxsec]')
end

% calculate number of events to add
% ---------------------------------
bg = 0;        % beginning of data
en = EEG.xmax; % end of data in sec
nu = floor(EEG.xmax/g.recurrence); % number of type 'X' events to add and epoch on

% bg = EEG.event(1).latency/EEG.srate;   % time in sec of first event
% en = EEG.event(end).latency/EEG.srate; % time in sec of last event
% nu = length((bg+g.recurrence):g.recurrence:(en-g.recurrence)); % number of 'X' events, one every 'g.recurrence' sec

if nu < 1
  error('specified recurrence_interval too long')
end

% print info on commandline
% -------------------------
eplength = g.limits(2)-g.limits(1);
fprintf('The input dataset will be split into %d epochs of %g s\n',nu,eplength);
fprintf('Epochs will overlap by %2.0f%%.\n',(eplength-g.recurrence)/eplength*100);

% insert events and urevents at the end of the current (ur)event tables
% ---------------------------------------------------------------------
fprintf('Inserting %d type ''%s'' events: ', nu, g.eventtype);
nevents = length(EEG.event);
% convert all event types to strings
for k = 1:nevents
    EEG.event(k).type = num2str(EEG.event(k).type);
end

nurevents = length(EEG.urevent);
for k = 1:nu
   if rem(k,40)
      fprintf('.')
   else
      fprintf('%d',k)
   end
   if k==40 | ( k>40 & ~rem(k-40,70))
     fprintf('\n');
   end

   EEG.event(nevents+k).type = g.eventtype;
   EEG.event(nevents+k).latency = g.recurrence*(k-1)*EEG.srate+1;

   EEG.urevent(nurevents+k).type = g.eventtype;
   EEG.urevent(nurevents+k).latency = g.recurrence*(k-1)*EEG.srate+1;
   EEG.event(nevents+k).urevent = nurevents+k;
end
fprintf('\n');

% sort the events based on their latency
% --------------------------------------
fprintf('Sorting the event table.\n');
EEG = pop_editeventvals( EEG, 'sort', {'latency' 0}); 

% split the dataset into epochs
% ------------------------------
if strcmpi(g.extractepochs, 'on')
    fprintf('Splitting the data into %d %2.2f-s epochs\n',nu,eplength); 
    setname = sprintf('%s - %g-s epochs', EEG.setname, g.recurrence);
    EEG = pop_epoch( EEG, { g.eventtype }, g.limits, 'newname', ...
                                      setname, 'epochinfo', 'yes');
    % baseline zero the epochs
    % ------------------------
    if ~isnan(g.rmbase)
        fprintf('Removing the pre %3.2f second baseline mean of each epoch.\n', g.rmbase);
        EEG = pop_rmbase( EEG, [g.limits(1) g.rmbase]);
    end
end;                              

