% eeg_interp() - interpolate data channels
%
% Usage: EEGOUT = eeg_interp(EEG, badchans, method);
%
% Inputs: 
%     EEG      - EEGLAB dataset
%     badchans - [integer array] indices of channels to interpolate.
%                For instance, these channels might be bad.
%                [chanlocs structure] channel location structure containing
%                either locations of channels to interpolate or a full
%                channel structure (missing channels in the current 
%                dataset are interpolated).
%     method   - [string] griddata method used for interpolation
%                (default is 'invdist')
%
% Output: 
%     EEGOUT   - data set with bad electrode data replaced by
%                interpolated data
%
% Author: Arnaud Delorme, CERCO, CNRS, Mai 2006

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
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

% $Log: not supported by cvs2svn $

function EEG = eeg_interp(ORIEEG, bad_elec, method)

    EEG = ORIEEG;
    if nargin < 2
        help eeg_interp;
        return;
    end;
    
    if nargin < 3
        method = 'invdist';
    end;

    if isstruct(bad_elec)
        
        % find missing channels
        % ---------------------
        if length(bad_elec) < length(EEG.chanlocs)
            bad_elec = [ EEG.chanlocs bad_elec ];
        end;
        if length(EEG.chanlocs) == length(bad_elec), return; end;
        
        lab1 = { bad_elec.labels };
        lab2 = { EEG.chanlocs.labels };
        [tmp badchans] = setdiff( lab1, lab2);
        fprintf('Found %d channels to interpolate\n');
        goodchans      = setdiff(1:length(bad_elec), badchans);
       
        % re-order good channels
        % ----------------------
        [tmp tmp2 neworder] = intersect( lab1, lab2 );
        [tmp2 ordertmp2] = sort(tmp2);
        neworder = neworder(ordertmp2);
        EEG.data = EEG.data(neworder, :, :);
        EEG.chanlocs = EEG.chanlocs(neworder); % not necessary
        if ~isempty(EEG.icasphere)
            EEG.icasphere = EEG.icasphere(:,neworder);
            EEG.icawinv   = pinv(EEG.icaweights*EEG.icasphere);
        end;
        
        % update EEG dataset (add blank channels)
        % ---------------------------------------
        if ~isempty(EEG.icasphere)
            if isempty(EEG.icachansind) || (length(EEG.icachansind) == EEG.nbchan)
                EEG.icachansind = goodchans; % this suppose that this is empty
            else
                error('Function not supported: cannot recompute ICA channel indices'); % just has to be programmed
            end;
        end;
        EEG.chanlocs             = bad_elec;
        tmpdata                  = zeros(length(bad_elec), size(EEG.data,2), size(EEG.data,3));
        tmpdata(goodchans, :, :) = EEG.data;
        EEG.data = tmpdata;
        EEG.nbchan = length(EEG.chanlocs);

    else
        badchans  = bad_elec;
        goodchans = setdiff(1:EEG.nbchan, badchans);
    end;

    % find non-empty good channels
    % ----------------------------
    nonemptychans = find(~cellfun('isempty', { EEG.chanlocs.theta }));
    [tmp indgood ] = intersect(goodchans, nonemptychans);
    goodchans = goodchans( sort(indgood) );
    
    % get theta, rad of electrodes
    % ----------------------------
    [xbad ,ybad]  = pol2cart([EEG.chanlocs( badchans).theta],[EEG.chanlocs( badchans).radius]);
    [xgood,ygood] = pol2cart([EEG.chanlocs(goodchans).theta],[EEG.chanlocs(goodchans).radius]);

    % scan data points
    % ----------------
    fprintf('Points:');
    for t=1:(size(EEG.data,2)*size(EEG.data,3)) % scan data points
        if mod(t,100) == 0, fprintf('%d ', t); end;
        if mod(t,1000) == 0, fprintf('\n'); end;
        %for c = 1:length(badchans)
        %   [h EEG.data(badchans(c),t)]= topoplot(EEG.data(goodchans,t),EEG.chanlocs(goodchans),'noplot', ...
        %        [EEG.chanlocs( badchans(c)).radius EEG.chanlocs( badchans(c)).theta]);
        %end;
        [Xi,Yi,EEG.data(badchans,t)] = griddata(ygood, xgood , EEG.data(goodchans,t)',...
                                                ybad, xbad, method); % interpolate data                                            
    end
    fprintf('\n');

    EEG = eeg_checkset(EEG);
