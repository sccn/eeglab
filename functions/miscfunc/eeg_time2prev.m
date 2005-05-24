% eeg_time2prev() - returns a vector giving, for each event of specified ("target") type(s), 
%                   the delay (in ms) since the preceding event (if any) of specified 
%                   ("previous") type(s). Requires the EEG.urevent structure, plus 
%                   EEG.event().urevent pointers to it. 
%
%           NOW SUPERCEDED BY eeg_context()
% Usage:
%           >> [delays,targets,urtargs,urprevs] = eeg_time2prev(EEG,{target},{previous});
% Inputs:
%  EEG      - structure containing an EEGLAB dataset
% {target}  - cell array of strings naming event type(s) of the specified target events 
%{previous} - cell array of strings naming event type(s) of the specified previous events 
%
% Outputs:
%  delays    - vector giving, for each event of a type listed in "target", the delay (in ms) 
%              since the last preceding event of a type listed in "previous" (0 if none such).
%  targets   - vector of indices of the "target" events in the event structure
%  urtargs   - vector of indices of the "target" events in the urevent structure
%  urprevs   - vector of indices of the "previous" events in the urevent structure (0 if none).
%
% Example:
% >> target = {'novel'};           % target event type 'novel'
% >> previous = {'novel', 'rare'}; % either 'novel' or 'rare'
% >> [delays,targets] = eeg_time2prev(EEG,target,previous);
%% Vector 'delays' now contains delays (in ms) from each 'novel' event to the previous
%% 'rare' OR 'novel' event, else 0 if none such. Vector 'targets' now contains the
%% 'novel' event indices.
%%
% Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, August 28, 2003

function [delays,targets,urtargets,urprevs] = eeg_time2prev(EEG,target,previous);

verbose = 1; % FLAG (1=on|0=off)
nevents = length(EEG.event);

%
%%%%%%%%%%%%%%%%%%%% Test input arguments %%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~iscell(target)
  error('2nd argument "target" must be a {cell array} of event type strings.');
  return
end
for k=1:length(target)
   if ~ischar([ target{k}{:} ])
       error('2nd argument "target" must be a {cell array} of event type strings.');
   end
end
if ~iscell(previous)
  error('3rd argument "previous" must be a {cell array} of event types.');
  return
end
for k=1:length(target) % make all target types strings
  if ~ischar(target{k})
    target{k} = num2str(target{k});
  end
end
for k=1:length(previous) % make all previous types strings
  if ~ischar(previous{k})
    previous{k} = num2str(previous{k});
  end
end
if ~isfield(EEG,'urevent')
  error('requires the urevent structure be present.');
  return
end
if length(EEG.urevent) < nevents
  error('WARNING: In dataset, number of urevents < number of events!?');
  return
end

% 
%%%%%%%%%%%%%%%%%%%% Initialize output arrays %%%%%%%%%%%%%%%%%%%%%%
%
delays    = zeros(1,nevents); % holds output times in ms
targets   = zeros(1,nevents); % holds indxes of targets
urtargets = zeros(1,nevents); % holds indxes of targets
urprevs   = zeros(1,nevents); % holds indxes of prevs

targetcount = 0; % index of current target
% Below:
%   idx   = current event index
%   uridx = current urevent index
%   tidx  = target type index
%   uidx  = previous urevent index
%   pidx  = previous type index
%
%%%%%%%%%%%for each event in the dataset %%%%%%%%%%%%%%%%%%%%
%
for idx = 1:nevents                                % for each event in the dataset
 %
 %%%%%%%%%%%%%%%%%%%%%%%% find target events %%%%%%%%%%%%%%%%%
 %
 uridx = EEG.event(idx).urevent;                   % find its urevent index
 istarget = 0;                                     % initialize target flag
 tidx = 1;                                         % initialize target type index
 while ~istarget & tidx<=length(target)            % for each potential target type
    if strcmpi(num2str(EEG.urevent(uridx).type),target(tidx))
      istarget=1;                                  % flag event as target
      targetcount = targetcount+1;                 % increment target count
      targets(targetcount) = idx;                  % save event index
      urtargets(targetcount) = uridx;              % save urevent index
      break                                        % stop target type checking
    else
      tidx = tidx+1;                               % else try next target type
    end
 end 

 if istarget  % if current event is a target type
  %
  %%%%%%%%%%%%%%%%%%% Find a 'previous' event %%%%%%%%%%%%%%%%%%
  %
  isprevious = 0;                                  % initialize previous flag
  uidx = uridx-1;                                  % begin with the previous urevent
  while uridx > 0
   pidx = 1;                                       % initialize previous type index
   while ~isprevious & pidx<=length(previous)      % for each previous event type
     if strcmpi(num2str(EEG.urevent(uidx).type),previous(pidx)) 
       isprevious=1;                               % flag 'previous' event
       urprevs(targetcount) = uidx;                % mark event as previous
       break                                       % stop previous type checking
     else
       pidx = pidx+1;                              % try next 'previous' event type
     end
   end  % pidx
   if isprevious
     break                                         % stop previous event checking
   else
     uidx = uidx-1;                                % keep checking for a 'previous' type event
   end % isprevious
  end % uidx
  %
  %%% Compute delay from current target to previous event %%%%%
  %
  if isprevious                                    % if type 'previous' event found
    delays(targetcount) = 1000/EEG.srate * ...
     (-1)*(EEG.urevent(uridx).latency - EEG.urevent(urprevs(targetcount)).latency);
  else
    delays(targetcount) = 0;                        % mark no previous event with 0
  end
  if verbose
    fprintf('%4d. (targ %s) %4d - (prev %s) %4d = %4.1f ms\n',...
                targetcount, targets(tidx),idx, ...
                             previous(pidx),urprevs(targetcount),...
                                          delays(targetcount));
  end % verbose
 end % istarget
end % event loop

%
%%%%%%%%%%%%%% Truncate the output arrays %%%%%%%%%%%%%%%%%%%%%%
%
if targetcount > 0
   targets   = targets(1:targetcount);
   urtargets = urtargets(1:targetcount);
   urprevs   = urprevs(1:targetcount);
   delays    = delays(1:targetcount);
else
   if verbose
     fprintf('eeg_time2prev(): No target type events found.\n')
   end
   targets   = [];
   urtargets = [];
   urprevs   = [];
   delays    = [];
end 
