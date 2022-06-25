% pop_mergeset() - Merge two or more datasets. If only one argument is given,
%                  a window pops up to ask for more arguments.
% Usage:
%   >> OUTEEG = pop_mergeset( ALLEEG ); % use a pop-up window
%   >> OUTEEG = pop_mergeset( ALLEEG, indices, keepall);
%   >> OUTEEG = pop_mergeset( INEEG1, INEEG2, keepall);
%
% Inputs:
%  INEEG1  - first input dataset
%  INEEG2  - second input dataset
%
% else
%  ALLEEG  - array of EEG dataset structures
%  indices - indices of EEG datasets to merge
%
%  keepall - [0|1] 0 -> remove, or 1 -> preserve, ICA activations
%            of the first dataset and recompute the activations
%            of the merged data {default: 0}
%
% Outputs:
%  OUTEEG  - merged dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% 01-25-02 reformated help & license -ad
% 01-26-02 change format for events and trial conditions -ad

function [INEEG1, com] = pop_mergeset( INEEG1, INEEG2, keepall);

com = '';
if nargin < 1
    help pop_mergeset;
    return;
end
if isempty(INEEG1)
    error('needs at least two datasets');
end
if nargin < 2 && length(INEEG1) == 1
    error('needs at least two datasets');
end

if nargin == 1
    uilist    = { { 'style' 'text' 'string' 'Dataset indices to merge' } ...
        { 'style' 'edit' 'string' '1' } ...
        { 'style' 'text' 'string' 'Preserve ICA weights of the first dataset ?' } ...
        { 'style' 'checkbox' 'string' '' } };
    res = inputgui( 'uilist', uilist, 'geometry', { [3 1] [3 1] }, 'helpcom', 'pophelp(''pop_mergeset'')');

    if isempty(res) return; end

    INEEG2  = eval( [ '[' res{1} ']' ] );
    keepall = res{2};
else
    if nargin < 3
        keepall = 0; % default
    end
end

fprintf('Merging datasets...\n');
if ~isstruct(INEEG2) % if INEEG2 is a vector of ALLEEG indices
    indices = INEEG2;

    NEWEEG = eeg_retrieve(INEEG1, indices(1)); % why abandoned?

    for index = 2:length(indices)
        INEEG2 = eeg_retrieve(INEEG1, indices(index));
        NEWEEG = pop_mergeset(NEWEEG, INEEG2, keepall); % recursive call
    end
    INEEG1 = NEWEEG;

else % INEEG is an EEG struct
    % check consistency
    % -----------------
    if INEEG1.nbchan ~= INEEG2.nbchan
        error('The two datasets must have the same number of channels');
    end
    if INEEG1.srate ~= INEEG2.srate
        error('The two datasets must have the same sampling rate');
    end
    if INEEG1.trials > 1 || INEEG2.trials > 1
        if INEEG1.pnts ~= INEEG2.pnts
            error('The two epoched datasets must have the same number of points');
        end
        if INEEG1.xmin ~= INEEG2.xmin
            INEEG2.xmin = INEEG1.xmin;
            fprintf('Warning: the two epoched datasets do not have the same time onset, adjusted');
        end
        if INEEG1.xmax ~= INEEG2.xmax
            INEEG2.xmax = INEEG1.xmax;
            fprintf('Warning: the two epoched datasets do not have the same time offset, adjusted');
        end
        
    end

    % repopulate epoch field if necessary
    % -----------------------------------
    if INEEG1.trials > 1 && INEEG2.trials == 1
        for iEvent = 1:length(INEEG2.event)
            INEEG2.event(iEvent).epoch = 1;
        end
    end
    if INEEG1.trials == 1 && INEEG2.trials > 1
        for iEvent = 1:length(INEEG1.event)
            INEEG1.event(iEvent).epoch = 1;
        end
    end
    
    % type of boundary event
    boundaryType = eeg_boundarytype(INEEG1, INEEG2);

    % Merge the epoch field
    % ---------------------
    if INEEG1.trials > 1 || INEEG2.trials > 1
        INEEGX = {INEEG1,INEEG2};
        for n = 1:2
            % make sure that both have an (appropriately-sized) epoch field
            % -------------------------------------------------------------
            if ~isfield(INEEGX{n},'epoch')
                INEEGX{n}.epoch = repmat(struct(),[1,INEEGX{n}.trials]);
            end
            % make sure that the epoch number is correct in each dataset
            % ----------------------------------------------------------
            if ~isempty(INEEGX{n}.epoch) && length(INEEGX{n}.epoch) ~= INEEGX{n}.trials
                disp('Warning: The number of trials does not match the length of the EEG.epoch field in one of');
                disp('         the datasets. Its epoch info will be reset and derived from the respective events.');
                INEEGX{n}.epoch = repmat(struct(),[1,INEEGX{n}.trials]);
            end
        end
        for n=1:2
            % purge all event-related epoch fields from each dataset (EEG.epoch.event* fields)
            % --------------------------------------------------------------------------------
            if isstruct(INEEGX{n}.epoch)
                fn = fieldnames(INEEGX{n}.epoch);
                INEEGX{n}.epoch = rmfield(INEEGX{n}.epoch,{fn{strmatch('event',fn)}});
                % copy remaining field names to the other dataset
                % -----------------------------------------------
                for f = fieldnames(INEEGX{n}.epoch)'
                    if ~isfield(INEEGX{3-n}.epoch,f{1})
                        INEEGX{3-n}.epoch(1).(f{1}) = [];
                    end
                end
                % after this, both sets have an epoch field with the appropriate number of items
                % and possibly some user-defined fields, but no event* fields.
            end
        end

        % concatenate epochs
        % ------------------
        if isstruct(INEEGX{1}.epoch) && isstruct(INEEGX{2}.epoch)
            if length(fieldnames(INEEGX{2}.epoch)) > 0
                INEEGX{1}.epoch(end+1:end+INEEGX{2}.trials) = orderfields(INEEGX{2}.epoch,INEEGX{1}.epoch);
            else
                INEEGX{1}.epoch(end+1:end+INEEGX{2}.trials) = INEEGX{2}.epoch;
            end
        end
        % and write back
        INEEG1 = INEEGX{1};
        INEEG2 = INEEGX{2};
        INEEGX = {};
    end
    
    % Concatenate data
    % ----------------
    if INEEG1.trials > 1 || INEEG2.trials > 1
        INEEG1.data(:,:,end+1:end+size(INEEG2.data,3)) = INEEG2.data(:,:,:);
    else
        INEEG1.data(:,end+1:end+size(INEEG2.data,2)) = INEEG2.data(:,:);
    end

    INEEG1.setname = 'Merged datasets';
    INEEG1trials = INEEG1.trials;
    INEEG2trials = INEEG2.trials;
    INEEG1pnts   = INEEG1.pnts;
    INEEG2pnts   = INEEG2.pnts;

    if INEEG1.trials > 1 || INEEG2.trials > 1 % epoched data
        INEEG1.trials  =  INEEG1.trials + INEEG2.trials;

    else % continuous data
        INEEG1.pnts = INEEG1.pnts + INEEG2.pnts;
    end

    if isfield(INEEG1, 'reject')
        INEEG1 = rmfield(INEEG1, 'reject' );
    end
    INEEG1.specicaact = [];
    INEEG1.specdata = [];

    if keepall == 0
        INEEG1.icaact = [];
        INEEG1.icawinv = [];
        INEEG1.icasphere = [];
        INEEG1.icaweights = [];
        if isfield(INEEG1, 'stats')
            INEEG1 = rmfield(INEEG1, 'stats' );
        end
    else
        INEEG1.icaact = [];
    end

    % concatenate events
    % ------------------
    if isempty(INEEG2.event) && INEEG2.trials == 1 && INEEG1.trials == 1 

        % boundary event
        % -------------
        disp('Inserting boundary event...');
        INEEG1.event(end+1).type    = boundaryType;     % add boundary event between datasets
        INEEG1.event(end  ).latency = INEEG1pnts+0.5; % make boundary halfway between last,first pts

        % check urevents
        % --------------
        if ~isfield(INEEG1, 'urevent'),
            INEEG1.urevent = [];
            fprintf('Warning: first dataset has no urevent structure.\n');
        end

        % add boundary urevent
        % --------------------
        disp('Inserting boundary urevent...');
        INEEG1.urevent(end+1).type    = boundaryType;

        if length(INEEG1.urevent) > 1 % if previous INEEG1 urevents
            INEEG1.urevent(end  ).latency = max(INEEG1pnts, INEEG1.urevent(end-1).latency)+0.5;
        else
            INEEG1.urevent(end  ).latency = INEEG1pnts+0.5;
        end

    else % is ~isempty(INEEG2.event)

        % concatenate urevents
        % --------------------
        if isfield(INEEG2, 'urevent')
            if ~isempty(INEEG2.urevent) && isfield(INEEG1.urevent, 'latency')

                % insert boundary event
                % ---------------------
                disp('Inserting boundary event...');
                INEEG1.urevent(end+1).type    = boundaryType;
                try
                    INEEG1.urevent(end  ).latency = max(INEEG1pnts, INEEG1.urevent(end-1).latency)+0.5;
                catch
                    % cko: sometimes INEEG1 has no events / urevents
                    INEEG1.urevent(end  ).latency = INEEG1pnts+0.5;
                end
                    

                % update urevent indices for second dataset
                % -----------------------------------------
                disp('Concatenating urevents...');
                orilen    = length(INEEG1.urevent);
                newlen    = length(INEEG2.urevent);
                INEEG2event = INEEG2.event;
                % update urevent index in INEEG2.event
                tmpevents = INEEG2.event;
                if ~isfield(tmpevents, 'urevent')
                    tmpevents(1).urevent = [];
                end
                nonemptymask = ~cellfun('isempty',{tmpevents.urevent});
                [tmpevents(nonemptymask).urevent] = celldeal(num2cell([INEEG2event.urevent]+orilen));
                INEEG2.event = tmpevents;
                % reserve space and append INEEG2.urevent
                INEEG1.urevent(orilen+newlen).latency = [];
                INEEG2urevent = INEEG2.urevent;
                tmpevents = INEEG1.urevent;
                for f = fieldnames(INEEG2urevent)'
                    [tmpevents((orilen+1):(orilen+newlen)).(f{1})] = INEEG2urevent.(f{1}); 
                end
                INEEG1.urevent = tmpevents;
            else
                INEEG1.urevent = [];
                INEEG2.urevent = [];
                fprintf('Warning: second dataset has empty urevent structure.\n');
            end
        end

        % concatenate events
        % ------------------
        disp('Concatenating events...');
        orilen = length(INEEG1.event);
        newlen = length(INEEG2.event);
        %allfields = fieldnames(INEEG2.event);

        % add discontinuity event if continuous
        % -------------------------------------
        if INEEG1trials  == 1 && INEEG2trials == 1
            disp('Adding boundary event...');
            INEEG1.event(end+1).type    = boundaryType; 
            INEEG1.event(end  ).latency = INEEG1pnts+0.5; 
            INEEG1.event(end  ).duration = NaN; 
            orilen = orilen+1;
%            eeg_insertbound(INEEG1.event, INEEG1.pnts, INEEG1pnts+1, 0); % +1 since 0.5 is subtracted
        end
        
        % ensure similar event structures
        % -------------------------------
        if ~isempty(INEEG2.event)
            if isstruct(INEEG1.event)
                for f = fieldnames(INEEG1.event)'
                    if ~isfield(INEEG2.event,f{1})
                        INEEG2.event(1).(f{1}) = [];
                    end
                end
            end
            if isstruct(INEEG2.event)
                for f = fieldnames(INEEG2.event)'
                    if ~isfield(INEEG1.event,f{1})
                        INEEG1.event(1).(f{1}) = [];
                    end
                end
            end
            INEEG2.event = orderfields(INEEG2.event, INEEG1.event);
        end

        % append
        % ------
        INEEG1.event(orilen + (1:newlen)) = INEEG2.event;
        INEEG2event = INEEG2.event;
        if isfield(INEEG1.event,'latency') && isfield(INEEG2.event,'latency')
            % update latency
            tmpevents = INEEG1.event;
            [tmpevents(orilen + (1:newlen)).latency] = celldeal(num2cell([INEEG2event.latency] + INEEG1pnts*INEEG1trials));
            INEEG1.event = tmpevents;
        end
        if isfield(INEEG1.event,'epoch') && isfield(INEEG2.event,'epoch')
            % update epoch index
            tmpevents = INEEG1.event;
            [tmpevents(orilen + (1:newlen)).epoch] = celldeal(num2cell([INEEG2event.epoch]+INEEG1trials));
            INEEG1.event = tmpevents;
        end

    end

    INEEG1.pnts = size(INEEG1.data,2);

    if ~isfield(INEEG1.event,'epoch') && ~isempty(INEEG1.event) && (size(INEEG1.data,3)>1 || ~isempty(INEEG1.epoch))
        INEEG1.event(1).epoch = [];
    end
    
    % rebuild event-related epoch fields
    % ----------------------------------
    disp('Reconstituting epoch information...');
    INEEG1.epoch = [];
    INEEG1 = eeg_checkset(INEEG1, 'eventconsistency');
end

% build the command
% -----------------
if nargout > 1
    if exist('indices') == 1
        com = sprintf('EEG = pop_mergeset( ALLEEG, [%s], %d);', int2str(indices), keepall);
    else
        com = sprintf('EEG = pop_mergeset( ALLEEG, EEG, %d);', keepall);
    end
 end

return


function varargout = celldeal(X)
varargout = X;
