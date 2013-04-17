% eeg_context() - returns (in output 'delays') a matrix giving, for each event of specified 
%                 ("target") type(s), the latency (in ms) to the Nth preceding and/or following
%                 urevents (if any) of specified ("neighbor") type(s). Return the target event 
%                 and urevent numbers, the neighbor urevent numbers, and the values of specified 
%                 urevent field(s) for each of the neighbor urevents. Uses the EEG.urevent 
%                 structure, plus EEG.event().urevent pointers to it. If epoched data, also 
%                 uses the EEG.epoch structure. For use in event-handling scripts and functions.
% Usage:
%             >>  [targs,urnbrs,urnbrtypes,delays,tfields,urnfields] = ...
%                          eeg_context(EEG,{targets},{neighbors},[positions],{fields},alltargs);
% Required input:
% EEG         - EEGLAB dataset structure containing EEG.event and EEG.urevent sub-structures
%
% Optional inputs:
% targets     - string or cell array of strings naming event type(s) of the specified target 
%               events {default | []: all events}
% neighbors   - string or cell array of strings naming event type(s) of the specified 
%               neighboring urevents {default | []: any neighboring events}.
% [positions] - int vector giving the relative positions of 'neighbor' type urevents to return. 
%               Ex: [-3 -2 -1 0 1 2 3] -> return the previous 3, current, and succeeding 3 
%               urevents of the specified {neighbor} types. [positions] values are arranged 
%               in ascending order before processing.  {default | []: 1 = first succeeding}
% fields      - string or cell array of strings naming one or more (ur)event field(s) to return 
%               values for neighbor urevents. {default: no field info returned}
% alltargs    - string ('all'|[]) if 'all', return information about all target urevents,
%               even those on which no epoch in the current dataset is centered. 
%               {default: [] -> only return information on epoch-centered target events} 
% Outputs:
%  targs      - size(ntargets,4) matrix giving the indices of target events in the event 
%               structure in column 1 and in the urevent structure in column 2. Column 3 gives 
%               the epoch number in which the target has latency 0 (else NaN if no such epoch).
%               The fourth column gives the index of the target type in the {targets} cell array.
%  urnbrs     - matrix of indices of "neighbor" events in the URevent structure (NaN if none).
%  urnbrtypes - int array giving the urnbrs event type indices in the {neighbor} cell array, 
%               else NaN if no such neighbor.  Ex: If nbr = {'square','rt'} (see below),
%               then urnbrtypes outputs are [1|2]; if nbr = {'rt'}, returns [1]s.
%  delays     - matrix giving, for each {targets} type event, the latency of the delay (in ms) 
%               from each target event to its neighbor urevents. Else, returns NaN when no 
%               neighbor event. Output matrix size: (ntargets,length(positions)). 
%  tfields    - real or cell array of values of the requested (ur)event field(s) for the target 
%               events. Values are the same type as the field values, else NaN if no such event.
%  urnfields  - real or cell array of values of the requested (ur)event field(s) for the neighbor 
%               urevents. Values are the same type as the field values, else NaN if no such event.
%               If > 1 field specified, a 3-D array or cell array (nevents,nnbrpos,nfields).
% Example:
%
% >> target   = 'square';                % target events are type 'square'
% >> nbr      = {'square','rt'};         % neighbor events are either 'square' or 'rt'
%
% >> [trgs,urnbrs,urnbrtypes,delays,tflds,urnflds] = ...
%                                          eeg_context(EEG,target,nbr,[-4 1],'position');
%    %
%    % Output 'delays' now contains latencies (in ms) from each 'square' target event to the 
%    % 4th preceding and 1st succeeding 'rt' OR 'square' urevent (else NaN when none such). 
%    % Outputs 'tfields' and 'urnflds' give the 'position' field values of target events and 
%    % neighbor urevents. Output 'urnbrtypes', the index of the type (1='square' or 2='rt') 
%    % of the ('urnbrs') neighbor urevents.
%    %
% >> trts      = find(trgs(:,4)==2);   % targets followed by an 'rt' (before any 'square')
% >> pos3      = find(trgfld = 3);     % targets with 'position'=3 (numeric field value).
% >> selevents = intersect_bc(trts,pos3); % target events by both criteria
% >> selepochs = trgs(selevents,3);    % epoch numbers centered on the selected target events
%
% Author: Scott Makeig, SCCN, Institute for Neural Computation, UCSD, March 27, 2004

% Edit History:
% 1/10/05 added 4th (type) field to output trgs; added to Example;
% 5/25/04 test for isnan(urevent.duration) to detect real break events -sm
% 3/27/04 made no-such-event return NaNs; added {target} and {neighbor} defaults -sm
% 3/28/04 added test for boundary urevents; renamed relidx as positions, lats as delays  -sm
% 3/29/04 reorganized output order, adding urnbrtypes -sm
% 3/29/04 changed urnbrtypes to int array -sm
% 3/31/04 ruled out searching 'boundary' events -sm
% 5/06/04 completed the function -sm
%
function [targs,ur_nbrs,ur_nbrtypes,delays,tfields,nfields] = eeg_context(EEG,targets,neighbors,positions,field,alltargs)

verbose     = 0;    % flag useful info printout (1=on|0=off)
debug_print = 0;    % flag overly verbose printout
breakwarning = 0;   % flag no pre-4.4 warning given

if nargin < 1
   help eeg_context
   return
end

if nargin< 6 | isempty(alltargs)
  alltargs = 0;
elseif strcmpi(alltargs,'all')
  alltargs = 1;
else
  error('alltargs argument must be ''all'' or [].')
end

if ~isstruct(EEG)
   error('first argument must be an EEG dataset structure');
end
if ~isfield(EEG,'event')
   error('No EEG.event structure found');
end
if ~isfield(EEG,'urevent')
   error('No EEG.urevent structure found');
end
if ~isfield(EEG.event,'urevent')
   error('No EEG.event().urevent field found');
end
      
if EEG.trials == 1 | ~isfield(EEG.event(1),'epoch')
  fprintf('Data are continuous: Returning info on all targets; no epoch info returned.\n')
  alltargs = 1;
  epochinfo = 0;
else
  epochinfo = 1;
end
if epochinfo & ~isfield(EEG,'epoch')
  error('No EEG.epoch information in this epoched dataset - run eeg_checkset()');
end
if epochinfo & ~isfield(EEG.epoch,'eventlatency')
  error('No EEG.epoch.eventlatency information in this epoched dataset');
end
if epochinfo & ~isfield(EEG.epoch,'event')
  error('No EEG.epoch.event information in this epoched dataset');
end

nevents   = length(EEG.event);
nurevents = length(EEG.urevent);
if length(EEG.urevent) < nevents
  fprintf('In this dataset there are more events than urevents. Check consistency.\n');
end
%
%%%%%%%%%%%%%%%%%% Substitute input defaults %%%%%%%%%%%%%%%%%%%%
%
if nargin < 5 | isempty(field)
  NO_FIELD = 1;        % flag no field variable output
end
if nargin  < 4 | isempty(positions)
  positions = 1;    % default: find next
end
if nargin < 3 | isempty(neighbors)
  neighbors = {'_ALL'};  % flag neighbors are all neighboring events
end
if nargin < 2 | isempty(targets)
  targets = {'_ALL'};  % flag targets are all events
end

%
%%%%%%%%%%%%% Test and adjust input arguments %%%%%%%%%%%%%%%%%%%%
%
if ischar(targets)
   targets = {targets};
end
if ~iscell(targets)
  error('2nd argument "targets" must be a {cell array} of event types.');
end
if ischar(neighbors)
   neighbors = {neighbors};
end
if ~iscell(neighbors)
  error('3rd argument "neighbors" must be a {cell array} of event types.');
end

for k=1:length(targets) % make all target types strings
  if ~ischar(targets{k})
    targets{k} = num2str(targets{k});
  end
end
for k=1:length(neighbors) % make all neighbor types strings
  if ~ischar(neighbors{k})
    neighbors{k} = num2str(neighbors{k});
  end
end

tmp = sort(positions);  % reorder positions in ascending order
if sum(tmp==positions) ~= length(positions)
    fprintf('eeg_context(): returning neighbors in ascending order: ');
    for k=1:length(tmp)
         fprintf('%d ',tmp(k));
    end
    fprintf('\n');
end
positions = tmp;

%
%%%%%%%%%%%%%% Prepare to find "neighbor" events %%%%%%%%%%%%%%%%%%
%
zeroidx = find(positions == 0);            % find 0's in positions vector
negidx = find(positions < 0);
negpos = positions(negidx);
if ~isempty(negpos)
  negpos = abs(negpos(end:-1:1));          % work backwards, make negpos positive
  negidx = negidx(end:-1:1);               % index into output ur_nbrs
end
nnegpos   = length(negpos);                % number of negative positions to search for

posidx = find(positions>0);                % work forwards
pospos = positions(posidx);
npospos = length(pospos);                  % number of positive positions to search for

% 
%%%%%%%%%%%%%%%%%%%% Initialize output arrays %%%%%%%%%%%%%%%%%%%%%%
%
npos      = length(positions);
delays    = NaN*zeros(nevents,npos); % holds inter-event intervals in ms 
                                     % else NaN when no neighboring event
targs     = NaN*zeros(nevents,1);    % holds indices of targets
targepochs= NaN*zeros(nevents,1);    % holds numbers of the epoch centered
                                     % on each target (or NaN if none such)
ur_trgs   = NaN*zeros(nevents,1);    % holds indices of targets
ur_nbrs   = NaN*zeros(nevents,npos); % holds indices of neighbors
ur_nbrtypes  = NaN*zeros(nevents,npos);  % holds {neighbors} type indices

cellfld = -1;  % flag no field output specified

if ~exist('NO_FIELD','var') % if field values asked for
  if ischar(field)
    if ~isfield(EEG.urevent,field)
     error('Specified field not found in urevent struct');
    end
    if ischar(EEG.urevent(1).(field)) ...
       | iscell(EEG.urevent(1).(field)) ...
          | isstruct(EEG.urevent(1).(field)), 
     tfields = cell(nevents,1);
     nfields = cell(nevents,npos);
     cellfld = 1;  % flag that field outputs are cell arrays
    else % is number
     tfields = NaN*zeros(nevents,1);
     nfields = NaN*zeros(nevents,npos);
     cellfld = 0;  % flag that field outputs are numeric arrays
    end 
    field = {field}; % make string field a cell array for uniformity
    nfieldnames = 1;
  elseif iscell(field)
    nfieldnames = length(field);
    for f = 1:nfieldnames
      if ~isfield(EEG.urevent,field{f})
       error('Specified field not found in urevent struct');
      end
      if ischar(EEG.urevent(1).(field{f})) ...
         | iscell(EEG.urevent(1).(field{f})) ...
           | isstruct(EEG.urevent(1).(field{f})),
        if f==1, 
           tfields = cell(nevents,nfieldnames);
           nfields = cell(nevents,npos,nfieldnames);
        end
        cellfld = 1;  % flag that field outputs are cell arrays
      else % is number
        if f==1
           tfields = NaN*zeros(nevents,nfieldnames);
           nfields = NaN*zeros(nevents,npos,nfieldnames);
        end
      end 
    end
    if cellfld == -1,
       cellfld = 0; % field value(s) must be numeric
    end
  else
     error('''field'' value must be string or cell array');
  end
end

targetcount = 0; % index of current target

% Below:
% evidx   = current event index
% uridx   = current urevent index
%  uidx   = neighbor urevent index
%  tidx   = target type index
%  nidx   = neighbor type index
%  pidx   = position index
%
%%%%%%%%%%%for each event in the dataset %%%%%%%%%%%%%%%%%%%%
%
wb=waitbar(0,'Computing event contexts...','createcancelbtn','delete(gcf)');
noepochs = 0;                                      % counter of targets that are not epoch-centered
targidx = zeros(nevents,1);

for evidx = 1:nevents  %%%%%% for each event in the dataset %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

waitbar(evidx/nevents);                            % update the waitbar fraction
if ~strcmp(EEG.event(evidx).type,'boundary')       % ignore boundary events (no urevents!)
 %
 %%%%%%%%%%%%%%%%%%%%%%%% find target events %%%%%%%%%%%%%%%%%
 %
 uridx = EEG.event(evidx).urevent;                 % find its urevent index
 if isempty(uridx)
   fprintf('eeg_context(): Data event %d does not point to an urevent.\n',evidx);
   delete(wb);
   return
 end
 istarget = 0;                                     % initialize target flag
 tidx = 1;                                         % initialize target type index
 %
 %%%%%%%%%%%%%%% cycle through target types %%%%%%%%%%%%%%%%%%
 %
 while ~istarget & tidx<=length(targets)           % for each potential target type
    uridxtype = EEG.urevent(uridx).type;
    if ~ischar(uridxtype), uridxtype = num2str(uridxtype); end
    if strcmpi(uridxtype,targets(tidx)) | strcmp(targets{1},'_ALL')         
                                                   % if is a target type
      istarget=1;                                  % flag event as target
      targetcount = targetcount+1;                 % increment target count
      %
      %%%%%%%%%%%%% find 0th neighbors (=targets themselves)
      %
      if ~isempty(zeroidx)  % if the target event is asked for in the nbrs array
        delays(targetcount,zeroidx) = 0;
        if ~exist('NO_FIELD','var')
            for f=1:nfieldnames
              if cellfld == 0
                nfields(targetcount,zeroidx,f) = EEG.urevent(uridx).(field{f});
              else % cellfld == 1
                nfields{targetcount,zeroidx,f} = EEG.urevent(uridx).(field{f});
              end
            end
        end
      end
      if epochinfo % if the data are epoched
            is0epoch = 0;
            for z = 1:length(EEG.event(evidx).epoch) % for each epoch the event is in
              ep = EEG.event(evidx).epoch(z);
              for e = 1:length(EEG.epoch(ep).event) % for each event in the epoch
                if EEG.epoch(ep).event(e) == evidx  % if current event
                    % js : added the if ~iscell loop
                  if ~iscell(EEG.epoch(ep).eventlatency) % i.e. not more than 1 eventtype in the epoch
                    if length(EEG.epoch(ep).eventlatency(e)) == 1
                      trglt = EEG.epoch(ep).eventlatency(e);    % get its epoch latency
                    else
                      trglt = EEG.epoch(ep).eventlatency(1); % this shouldn't happen
                      fprintf('EEG.epoch(%d).eventlatency(%d) length > 1 ??\n',ep,e);
                    end
                  else
                    if length(EEG.epoch(ep).eventlatency{e}) == 1  
                      trglt = EEG.epoch(ep).eventlatency{e};    % get its epoch latency 
                    else
                      trglt = EEG.epoch(ep).eventlatency{e}(1); % this shouldn't happen
                      fprintf('EEG.epoch(%d).eventlatency{%d} length > 1 ??\n',ep,e);
                    end;
                  end
                  if trglt == 0
                      targepochs(targetcount) = ep;
                      is0epoch = 1;
                      break;
                  end
                end
              end
              if is0epoch 
                break; 
              end
            end % for
            if ~is0epoch, noepochs = noepochs+1; end 
      end

       targs(targetcount) = evidx;                   % save event index
       ur_trgs(targetcount) = uridx;                 % save urevent index
       if ~exist('NO_FIELD','var')
         for f=1:nfieldnames
           if cellfld ==0
            tfields(targetcount,f) = EEG.urevent(uridx).(field{f});
           elseif cellfld == 1
            tfields{targetcount,f} = EEG.urevent(uridx).(field{f});
           end
         end
         break                                       % stop target type checking
       end
    else % next target type
       tidx = tidx+1;                                % else try next target type
    end % if is target
 end % while ~istarget

 if istarget  % if current event is a target type
  targidx(targetcount) = tidx;                       % save index of its type within targets
  if ~isempty(negpos)
   %
   %%%%%%%%%%%%%%%%%% find previous neighbor urevents %%%%%%%%%%%%
   %
   uidx = uridx-1;                                   % begin with the previous urevent
   npidx  = 1;                                       % index into negpos
   curpos = 1;                                       % current (negative) position
   seekpos = negpos(npidx);                          % begin with first negpos position
   while uidx > 0 & npidx <= nnegpos                 % search through previous urevents
     uidxtype = EEG.urevent(uidx).type;
     if ~ischar(uidxtype), uidxtype = num2str(uidxtype); end
     if strcmpi(uidxtype,'boundary')  % flag boundary urevents
        if ~isfield(EEG.urevent,'duration') ...
          | ( isnan(EEG.urevent(uidx).duration) ...
             | isempty(EEG.urevent(uidx).duration))  % pre-v4.4 or real break urevent 
                                                     %   (NaN duration)
             if ~isfield(EEG.urevent,'duration') ... % pre version-4.4 dataset
                    & breakwarning == 0
               fprintf('Pre-v4.4 boundary urevent found - duration field not defined.');
               breakwarning = 1;
             end
        end
        break           % don't search for neighbors across a boundary urevent
     end
     isneighbor = 0;                                 % initialize neighbor flag
     nidx = 1;                                       % initialize neighbor type index
     %
     %%%%%%%%%% cycle through neighbor types %%%%%%%%%%%%%
     %
     while ~isneighbor & nidx<=length(neighbors)     % for each neighbor event type
       if strcmpi(uidxtype,neighbors(nidx)) | strcmp(neighbors,'_ALL')
         isneighbor=1;                               % flag 'neighbors' event
         curpos = curpos+1;
         %
         %%%%%%%%%%%%%%% if an event in one of the specified positions %%%%%
         %
         if curpos-1 == seekpos
           delays(targetcount,negidx(npidx)) = 1000/EEG.srate * ...
                                (EEG.urevent(uidx).latency - EEG.urevent(uridx).latency);
                                % return negative latencies for negpos events
           ur_nbrs(targetcount,negidx(npidx))=uidx;  % mark urevent as neighbor
           ur_nbrtypes(targetcount,negidx(npidx)) = nidx;
           if ~exist('NO_FIELD','var')
             for f=1:nfieldnames
               if cellfld ==0
                if nfieldnames > 1
                   if ~isempty(EEG.urevent(uidx).(field{f}))
                     nfields(targetcount,negidx(npidx),f) = EEG.urevent(uidx).(field{f});
                   else
                     nfields(targetcount,negidx(npidx),f) = NaN;
                   end
                elseif ~isempty(EEG.urevent(uidx).(field{f}))
                   nfields(targetcount,negidx(npidx)) = EEG.urevent(uidx).(field{f});
                else
                   nfields(targetcount,negidx(npidx)) = NaN;
                end
               elseif cellfld == 1
                if nfieldnames > 1
                  nfields{targetcount,negidx(npidx),f} = EEG.urevent(uidx).(field{f});
                else
                  nfields{targetcount,negidx(npidx)} = EEG.urevent(uidx).(field{f});
                end
               end
             end
           end
           npidx = npidx+1;  % look for next negpos position
           if npidx<=nnegpos
              seekpos = negpos(npidx);               % new seek position
           end
         end % if seekpos
         break                                       % stop neighbors type checking
       else
         nidx = nidx+1;                              % try next 'neighbors' event type
       end
     end  % nidx - neighbor-type testing loop
     %
     %%%%%%%%%%%%%%% find preceding neighbor event %%%%%%%%%%%%%%
     %
     uidx = uidx-1;                     % keep checking for a 'neighbors' type event
   end % while uidx - urevent type checking
  end % if negpos 

  if ~isempty(pospos)
   %
   %%%%%%%%%%%%%%% find succeeding position urevents %%%%%%%%%%%%
   %
   uidx = uridx+1;                                   % begin with the succeeding urevent
   ppidx  = 1;                                       % index into pospos
   curpos = 1;                                       % current (positive) position
   seekpos = pospos(ppidx);                          % begin with first pospos position
   while uidx <= nurevents  & ppidx <= npospos       % search through succeeding urevents
     isneighbor = 0;                                 % initialize neighbor flag
     uidxtype = EEG.urevent(uidx).type;
     if ~ischar(uidxtype), uidxtype = num2str(uidxtype); end
     if strcmpi(uidxtype,'boundary')  % flag boundary events
        if ~isfield(EEG.urevent,'duration') ...
          | ( isnan(EEG.urevent(uidx).duration) ...
              | isempty(EEG.urevent(uidx).duration))  % pre-v4.4 or real break urevent 
                                                      %   (NaN duration)
             if ~isfield(EEG.urevent,'duration') ...  % pre version-4.4 dataset
                    & breakwarning == 0
               fprintf('Pre-v4.4 boundary urevent found - duration field not defined.');
               breakwarning = 1;
             end
        end
        break           % don't search for neighbors across a boundary urevent
     end
     pidx = 1;                                        % initialize neighbor type index
     %
     %%%%%%%%%% cycle through neighbor types %%%%%%%%%%%%%
     %
     while ~isneighbor & pidx<=length(neighbors)      % for each neighbor event type
       if strcmpi(uidxtype,neighbors(pidx)) ...
                          | strcmp(neighbors,'_ALL')
         isneighbor=1;                                % flag 'neighbors' event
         curpos = curpos+1;
         %
         %%%% if an event in one of the specified positions %%%%%
         %
         if curpos-1 == seekpos
           ur_nbrs(targetcount,posidx(ppidx))=uidx;   % mark urevent as neighbor
           % ur_nbrtypes{targetcount,posidx(ppidx)} = EEG.urevent(uidx).type; % note its type
           ur_nbrtypes(targetcount,posidx(ppidx)) = pidx; % note its type
           delays(targetcount,posidx(ppidx)) = 1000/EEG.srate * ...
                                (EEG.urevent(uidx).latency - EEG.urevent(uridx).latency);
                                % return positive latencies for pospos events
           if ~exist('NO_FIELD','var')
             for f=1:nfieldnames
               if cellfld ==0
                  if nfieldnames > 1
                    if ~isempty(EEG.urevent(uidx).(field{f}))
                      nfields(targetcount,posidx(ppidx),f) = EEG.urevent(uidx).(field{f});
                    else
                      nfields(targetcount,posidx(ppidx),f) = NaN;
                    end
                  elseif ~isempty(EEG.urevent(uidx).(field{f}))
                    nfields(targetcount,posidx(ppidx)) = EEG.urevent(uidx).(field{f});
                  else
                    nfields(targetcount,posidx(ppidx)) = NaN;
                  end
               else % if cellfld == 1
                  if nfieldnames > 1
                    nfields{targetcount,posidx(ppidx),f} = EEG.urevent(uidx).(field{f});
                  else
                    nfields{targetcount,posidx(ppidx)} = EEG.urevent(uidx).(field{f});
                  end
               end
             end
           end
           ppidx = ppidx+1;
           if ppidx<=npospos
              seekpos = pospos(ppidx);                % new seek position
           end
           break                                      % stop neighbors type checking
         end % if seekpos
         break
       else
         %
         %%%%%%%%%%%% find succeeding neighbor event %%%%%%%%%%%%
         %
         pidx = pidx+1;                               % try next 'neighbors' event-type
       end
     end  % pidx - neighbor-type testing loop
     uidx = uidx+1;                     % keep checking for 'neighbors' type urevents
   end % uidx - urevent type checking
  end % if pospos 
  %
  %%%%%%%%%%%%%%% debug mode info printout %%%%%%%%%%%%%%%%%%%%%%
  %
  if debug_print
    fprintf('%d. ',targetcount)
    if targetcount<1000, fprintf(' '); end
    if targetcount<100, fprintf(' '); end
    if targetcount<10, fprintf(' '); end;
    if uidx > 1
      %fprintf('event %-4d ttype %s - delays: ',evidx,num2str(EEG.urevent(evidx).type));
     for k=1:npos
      fprintf('(%d) ',ur_nbrs(targetcount,k));
      if ur_nbrs(targetcount,k)<1000, fprintf(' '); end
      if ur_nbrs(targetcount,k)<100, fprintf(' '); end
      if ur_nbrs(targetcount,k)<10, fprintf(' '); end;
      fprintf('%2.0f ',delays(targetcount,k));
     end
     if ~exist('NO_FIELD','var')
       if cellfld == 0 % numeric field values
        fprintf('fields: ')
        for f=1:nfieldnames
          fprintf('%-5g - ',tfields(targetcount,f))
          for k=1:npos
            fprintf('%-5g ',nfields(targetcount,k,f));
          end
        end
       elseif cellfld == 1  % cell array field values
        fprintf('fields: ')
        for f=1:nfieldnames
          if ischar(EEG.urevent(1).(field{f}))
           fprintf('%-5g -',tfields{targetcount,f})
           for k=1:npos
              fprintf('%-5g ',nfields{targetcount,k,f});
           end % k
          end % ischar
        end % f
       end % cellfield
     end % ~NO_FIELD
    end % uidx > 1
    fprintf('\n');
  end % debug_print
 end % istarget
 %
 %%%%%%%%%%%%%%%%% find next target event %%%%%%%%%%%%%%%%%%%%%%
 %
 end % if not 'boundary' event
 evidx = evidx+1;      % continue event checking
end % event loop

%
%%%%%%%%% delete watibar %%%%%%%%
%
if ishandle(wb), delete(wb); end;

if ~alltargs
  fprintf('Returning info on the %d of %d target events that have epochs centered on them.\n',...
                  targetcount-noepochs,targetcount);
else
  fprintf('Returning info on all %d target events (%d have epochs centered on them).\n',...
                  targetcount,targetcount-noepochs);
end

if debug_print
  fprintf('---------------------------------------------------\n');
  fprintf('ur#   event   #  ttype targtype - delays (urnbr) ms\n');
end
%
%%%%%%%%%%%%%% Truncate the output arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if targetcount > 0
   targs   = [targs(1:targetcount) ur_trgs(1:targetcount) ...
                     targepochs(1:targetcount) targidx(1:targetcount)]; % 4-column output
   ur_nbrs   = ur_nbrs(1:targetcount,:);
   delays    = delays(1:targetcount,:);

   epcenttargs = find(~isnan(targs(:,3))); % find targets that have an epoch centered on them
   if ~alltargs
      targs = targs(epcenttargs,:);
      ur_nbrs = ur_nbrs(epcenttargs,:);
      delays = delays(epcenttargs,:);
   end

   if ~exist('NO_FIELD','var')
     if cellfld == 0
       tfields = tfields(1:targetcount,:);
       nfields = nfields(1:targetcount,:,:);
     elseif cellfld == 1
       tfields = tfields(1:targetcount,:);
       nfields = nfields(1:targetcount,:,:);
     end
     if ~alltargs
        tfields = tfields(epcenttargs,:);
        nfields = nfields(epcenttargs,:,:);
     end
   else % NO_FIELD
     tfields = [];
     nfields = [];
   end
   ur_nbrtypes = ur_nbrtypes(1:targetcount,:);
   if ~alltargs
      ur_nbrtypes = ur_nbrtypes(epcenttargs,:);
   end
else % return nulls if no targets found
   if verbose
     fprintf('eeg_context(): No target type events found.\n')
   end
   delays  = [];
   targs   = [];
   ur_nbrs = [];
   tfields = [];
   nfields = [];
end 

