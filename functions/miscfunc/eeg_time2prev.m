% eeg_time2prev() - returns a vector giving, for each event of a specified ("target") type, 
%                   the delay (in ms) since the past preceding event (if any) of an event of 
%                   a specified ("previous") type. Requires the EEG.urevent structure, plus 
%                   EEG.event().urevent pointers to it. 
% Usage:
%      >> [delays,targets,urtargs,urprevs] = eeg_time2prev(EEG,{target},{previous});
% Inputs:
%  EEG      - structure containing an EEGLAB dataset
%  target   - cell array of strings naming event type(s) of the specified target events 
%  previous - cell array of strings naming event type(s) of the specified previous events 
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
%%
%% Vector delays now contains delays in ms from each 'novel' event to the previous
%% 'rare' or 'novel' event (else 0 when none such). Vector targets now contains the
%% 'novel' event indices.
%%
% Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, August 28, 2003

function [delays,targets,urtargets,urprevs] = eeg_time2prev(EEG,target,previous);

verbose = 1; % 1=on, 0=off
targetcount = 0;
nevents   = length(EEG.urevent);

if ~iscell(target)
  error('2nd arg "target" must be a cell array.\n');
  return
end
if ~iscell(previous)
  error('3rd arg "previous" must be a cell array.\n');
  return
end
if ~isfield(EEG,'urevent')
  error('requires the urevent structure be present.\n');
  return
end
if length(EEG.urevent) < length(EEG.event)
  error('number of urevents < number of events!?\n');
  return
end

delays    = zeros(1,nevents); % holds output times in ms
targets   = zeros(1,nevents); % holds indxes of targets
urtargets = zeros(1,nevents); % holds indxes of targets
urprevs   = zeros(1,nevents); % holds indxes of prevs

for idx = 1:length(EEG.event)
 uridx = EEG.event(idx).urevent;
 istarget = 0; tidx = 1;
 while istarget==0 & tidx<=length(target)
    if strcmpi(EEG.urevent(uridx).type,target(tidx))
      istarget=1;
      targetcount = targetcount+1;
      targets(targetcount) = idx;
      urtargets(targetcount) = uridx;
    end
    tidx = tidx+1;
 end 

 if istarget 
  if urprevs(targetcount) > 0   % fill delays with latency difference in ms
    delays(targetcount) = 1000/EEG.srate * ...
     (-1)*(EEG.urevent(uridx).latency - EEG.urevent(urprevs(targetcount)).latency);
  elseif targetcount>1 & urprevs(targetcount-1)>0
    urprevs(targetcount) = urprevs(targetcount-1);
    delays(targetcount) = 1000/EEG.srate * ...
     (-1)*(EEG.urevent(uridx).latency - EEG.urevent(urprevs(targetcount)).latency);
  end
  if verbose
    fprintf('%4d. target event %4d - previous event %4d = %4.1f ms\n',...
               targetcount,idx,urprevs(targetcount),delays(targetcount));
  end
 end

 isprevious = 0; pidx = 1;
 while ~isprevious & pidx<=length(previous)
    if strcmpi(EEG.urevent(uridx).type,previous(pidx))
      isprevious=1;
      urprevs(targetcount+1) = uridx;
    end
    pidx = pidx+1;
 end 
end

if targetcount > 0
   targets   = targets(1:targetcount);
   urtargets = urtargets(1:targetcount);
   urprevs   = urprevs(1:targetcount);
   delays    = delays(1:targetcount);
else
   fprintf('eeg_time2prev(): No target type events found.\n')
   targets   = [];
   urtargets = [];
   urprevs   = [];
   delays    = [];
end 
