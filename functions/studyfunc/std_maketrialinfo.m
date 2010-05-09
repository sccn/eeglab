function STUDY = std_maketrialinfo(STUDY, ALLEEG);

%% Make trial info
for index = 1:length(ALLEEG)
    eventlat = abs(eeg_point2lat( [ ALLEEG(index).event.latency ], [ ALLEEG(index).event.epoch ], ALLEEG(index).srate, [ALLEEG(index).xmin ALLEEG(index).xmax]));
    events   = ALLEEG(index).event;
    ff = fieldnames(events);
    ff = setdiff(ff, { 'latency' 'urevent' 'epoch' });
    trialinfo = [];
    
    % process time locking event fields
    % ---------------------------------
    indtle    = find(eventlat < 0.0005);
    epochs    = [ events(indtle).epoch ];
    if length(epochs) ~= ALLEEG(index).trials
        error('Not the same number of time-locking events as trials');
    end;
    commands = {};
    for f = 1:length(ff)
        eval( [ 'eventvals = {events(indtle).' ff{f} '};' ]);
        if isnumeric(eventvals{1})
            eventvals = cellfun(@num2str, eventvals, 'uniformoutput', false);
        end;
        commands = { commands{:} ff{f} eventvals };
    end;
    trialinfo = struct(commands{:});
    
%    % same as above but 10 times slower
%     for e = 1:length(ALLEEG(index).event)
%         if eventlat(e) < 0.0005 % time locking event only
%             epoch = events(e).epoch;
%             for f = 1:length(ff)
%                 fieldval  = getfield(events, {e}, ff{f});
%                 trialinfo = setfield(trialinfo, {epoch}, ff{f}, fieldval);
%             end;
%         end;
%     end;
    STUDY.datasetinfo(index).trialinfo = trialinfo;
end;

    