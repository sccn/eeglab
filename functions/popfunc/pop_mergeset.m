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

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.43  2007/04/27 22:26:10  arno
% trying to fix event merge
%
% Revision 1.42  2006/04/14 17:40:11  arno
% fixing INEEG1
%
% Revision 1.41  2006/04/08 04:02:54  toby
% corrected problems with events and epochs
%
% Revision 1.40  2006/04/08 01:58:55  toby
% removed some extraneous code used for bug checking
%
% Revision 1.39  2005/11/29 22:57:26  toby
% *** empty log message ***
%
% Revision 1.37  2005/11/10 23:42:03  toby
% edit to be compatible with Matlab 5
%
% Revision 1.35  2005/11/02 07:18:10  scott
% disp() statements
%
% Revision 1.34  2005/10/14 18:55:27  hilit
% fixed a bug that disabled merging epoched datasets
%
% Revision 1.33  2005/05/12 16:14:58  arno
% fix event latency boundary
%
% Revision 1.32  2004/12/08 18:29:20  arno
% creating event table if does not exist to insert boundary event
%
% Revision 1.31  2004/11/18 23:20:25  arno
% fixed urevent boundary latency
%
% Revision 1.30  2004/11/15 22:36:52  arno
% latency of urevent boundary
%
% Revision 1.29  2004/11/09 23:51:00  arno
% boundary in urevent
%
% Revision 1.28  2004/06/08 17:32:12  arno
% update new call to eeg_insertbound
%
% Revision 1.27  2004/05/14 22:14:54  arno
% same
%
% Revision 1.26  2004/05/14 22:10:18  arno
% new eeg_insertbound call
%
% Revision 1.25  2004/05/06 21:58:27  arno
% setfields of new events
%
% Revision 1.24  2004/05/06 17:47:22  arno
% debug add boundary event
%
% Revision 1.23  2004/05/06 17:42:14  arno
% add urevent boundary
%
% Revision 1.22  2004/03/19 19:15:58  arno
% merging data differently to save mem
%
% Revision 1.21  2004/03/17 23:10:16  arno
% fixing urevents
%
% Revision 1.20  2004/02/21 03:14:53  arno
% same
%
% Revision 1.19  2004/02/21 03:14:05  arno
% fixing urevent boundaries
%
% Revision 1.18  2004/02/21 00:23:04  arno
% urevent update
%
% Revision 1.17  2004/02/20 23:53:25  arno
% urevent merging
%
% Revision 1.16  2003/12/19 00:53:17  arno
% adding dicontinuity dataset
%
% Revision 1.15  2003/11/18 16:56:33  scott
% text
%
% Revision 1.14  2003/11/18 16:50:36  scott
% text
%
% Revision 1.13  2003/01/24 19:43:50  arno
% same
%
% Revision 1.12  2003/01/24 19:41:35  arno
% debugging ?
% /
%
% Revision 1.11  2003/01/16 18:30:35  arno
% allowing merging between dataset with different event structures
%
% Revision 1.10  2002/09/04 22:17:19  luca
% improving epoch merge checking -arno
%
% Revision 1.9  2002/08/13 23:58:35  arno
% update error message
%
% Revision 1.8  2002/08/12 02:37:13  arno
% inputdlg2
%
% Revision 1.7  2002/06/25 02:20:29  arno
% preserving epoch information
%
% Revision 1.6  2002/06/25 00:52:27  arno
% debuging ICA info copy
%
% Revision 1.5  2002/04/23 20:08:39  arno
% full reprogramming of the function for standalone
%
% Revision 1.4  2002/04/21 01:09:24  scott
% edited help msg -sm
%
% Revision 1.3  2002/04/18 20:03:45  arno
% retrIeve
%
% Revision 1.2  2002/04/10 21:32:26  arno
% debuging event concatenation
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 01-26-02 change format for events and trial conditions -ad

function [INEEG1, com] = pop_mergeset( INEEG1, INEEG2, keepall);

com = ''; 
if nargin < 1
	help pop_mergeset;
	return;
end;
if isempty(INEEG1)
	error('needs at least two datasets');
end;
if nargin < 2 & length(INEEG1) == 1
	error('needs at least two datasets');
end;

if nargin == 1
	promptstr    = { 'Dataset indices to merge', ...
					 'Preserve ICA weights of the first dataset ?' };
	inistr       = { '1', 'no' };
	result       = inputdlg2( promptstr, 'Merge datasets -- pop_mergeset()', 1,  inistr, 'pop_mergeset');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
   
    INEEG2  = eval( [ '[' result{1} ']' ] );
    switch lower(result{2})
    	case 'yes', keepall = 1;
    	otherwise, keepall = 0;
    end;
else
	if nargin < 3
		keepall = 0; % default
	end;	
end;

fprintf('Merging datasets...\n');		

if ~isstruct(INEEG2) % if INEEG2 is a vector of ALLEEG indices	
	indices = INEEG2;
	if length(indices) < 2
		error('needs at least two datasets');
	end;

	NEWEEG = eeg_retrieve(INEEG1, indices(1)); % why abandoned?

	for index = 2:length(indices)
		INEEG2 = eeg_retrieve(INEEG1, indices(index));
        NEWEEG = pop_mergeset(NEWEEG, INEEG2, keepall); % recursive call
	end;
	INEEG1 = NEWEEG;

else % INEEG is an EEG struct
	% check consistency
	% -----------------
	if INEEG1.nbchan ~= INEEG2.nbchan
		error('The two datasets must have the same number of channels');
	end;	
	if INEEG1.srate ~= INEEG2.srate
		error('The two datasets must have the same sampling rate');
	end;	
	if INEEG1.trials > 1 | INEEG2.trials > 1
		if INEEG1.pnts ~= INEEG2.pnts
			error('The two epoched datasets must have the same number of points');
		end;
		if INEEG1.xmin ~= INEEG2.xmin
                        INEEG2.xmin = INEEG1.xmin;
			fprintf('Warning: the two epoched datasets do not have the same time onset, adjusted');
		end;
		if INEEG1.xmax ~= INEEG2.xmax
                        INEEG2.xmax = INEEG1.xmax;
			fprintf('Warning: the two epoched datasets do not have the same time offset, adjusted');
		end;
	end;	

    % Concatenate data
    % ----------------
    if INEEG1.trials > 1 | INEEG2.trials > 1
        INEEG1.data(:,:,end+1:end+size(INEEG2.data,3)) = INEEG2.data;
    else
        INEEG1.data(:,end+1:end+size(INEEG2.data,2)) = INEEG2.data;
    end;

    INEEG1.setname = 'Merged datasets';
    INEEG1trials = INEEG1.trials;
    INEEG2trials = INEEG2.trials;
    INEEG1pnts   = INEEG1.pnts;
    INEEG2pnts   = INEEG2.pnts;

	if INEEG1.trials > 1 | INEEG2.trials > 1 % epoched data
		INEEG1.trials  =  INEEG1.trials + INEEG2.trials;

	else % continuous data
		INEEG1.pnts = INEEG1.pnts + INEEG2.pnts;
	end;
	
	if isfield(INEEG1, 'reject')
		INEEG1 = rmfield(INEEG1, 'reject' );
	end;
	INEEG1.specicaact = [];
	INEEG1.specdata = [];

	if keepall == 0
		INEEG1.icaact = [];
		INEEG1.icawinv = [];
		INEEG1.icasphere = [];
		INEEG1.icaweights = [];
		if isfield(INEEG1, 'stats')
			INEEG1 = rmfield(INEEG1, 'stats' );
		end;
	else
		INEEG1.icaact = [];
	end;

	% concatenate events
	% ------------------
	if isempty(INEEG2.event)

        % boundary event
        % -------------
        disp('Inserting boundary event...');
        INEEG1.event(end+1).type    = 'boundary';     % add boundary event between datasets
        INEEG1.event(end  ).latency = INEEG1pnts+0.5; % make boundary halfway between last,first pts
        
        % check urevents
        % --------------
        if ~isfield(INEEG1, 'urevent'), 
   			INEEG1.urevent = []; 
            fprintf('Warning: first dataset has no urevent structure.\n');
		end;

        % add boundary urevent
        % --------------------
        disp('Inserting boundary urevent...');
        INEEG1.urevent(end+1).type    = 'boundary';  

        if length(INEEG1.urevent) > 1 % if previous INEEG1 urevents
            INEEG1.urevent(end  ).latency = max(INEEG1pnts, INEEG1.urevent(end-1).latency)+0.5;
        else
            INEEG1.urevent(end  ).latency = INEEG1pnts+0.5;
        end;

    else % is ~isempty(INEEG2.event)

        % concatenate urevents
        % --------------------
        if isfield(INEEG2, 'urevent') 
		  if ~isempty(INEEG2.urevent) 
            
            % insert boundary event
            % ---------------------
            disp('Inserting boundary event...');
            INEEG1.urevent(end+1).type    = 'boundary';
            INEEG1.urevent(end  ).latency = max(INEEG1pnts, INEEG1.urevent(end-1).latency)+0.5;
            
            % update urevent indices for second dataset
            % -----------------------------------------
            disp('Concatenating urevents...');
            orilen    = length(INEEG1.urevent);
            for e=1:length(INEEG2.event)
                INEEG2.event(e).urevent = INEEG2.event(e).urevent + orilen;
            end;
            
            allfields = fieldnames(INEEG2.urevent);
            for i=1:length( allfields )
                for e=1:length(INEEG2.urevent)
                    tmpval = getfield(INEEG2.urevent, { e }, allfields{i});
                    INEEG1.urevent = setfield(INEEG1.urevent, {orilen + e}, allfields{i}, tmpval);
                end
            end
		  else
            fprintf('Warning: second dataset has empty urevent structure.\n');
          end 
        end;

        % concatenate events
        % ------------------
        disp('Concatenating events...');
        orilen = length(INEEG1.event);
        %allfields = fieldnames(INEEG2.event);

        % ensure similar event structures
        % -------------------------------
        if ~isempty(INEEG2.event)
            fields1 = lower(fieldnames(INEEG1.event));
            fields2 = lower(fieldnames(INEEG2.event));
            if length(fields1) > length(fields2)
                for index = 1:length(fields1)
                    if isempty(strmatch(fields1{index}, fields2))
                        INEEG1.event = setfield( INEEG1.event, { orilen + 1}, fields1{index}, []);
                    end;
                end;
            elseif length(fields1) < length(fields2)
                for index = 1:length(fields2)
                    if isempty(strmatch(fields2{index}, fields1))
                        INEEG2.event = setfield( INEEG2.event, { 1 }, fields2{index}, []);
                    end;
                end;
            end;
        end;
        
        for e=1:length(INEEG2.event)
            INEEG1.event(orilen + e) = INEEG2.event(e);
            if isfield(INEEG1.event,'latency') & isfield(INEEG2.event,'latency')
               INEEG1.event(orilen + e).latency = INEEG2.event(e).latency + INEEG1pnts * INEEG1trials;
            end
            if isfield(INEEG1.event,'epoch') & isfield(INEEG2.event,'epoch')
               INEEG1.event(orilen + e).epoch = INEEG2.event(e).epoch + INEEG1trials;
            end
         end

        INEEG1.epoch = []; % epoch info regenerated below by 'eventconsistency' in eeg_checkset()
  
        % add discontinuity event if continuous
        % -------------------------------------
        if INEEG1trials  == 1 & INEEG2trials == 1
            disp('Adding boundary event...');
            INEEG1.event = eeg_insertbound(INEEG1.event, INEEG1.pnts, INEEG1pnts+1, 0); % +1 since 0.5 is subtracted
        end;

	end;
        
	if isfield(INEEG1, 'epoch') && isfield(INEEG2, 'epoch') && ~isempty(INEEG1.epoch) && ~isempty(INEEG2.epoch)
		try 
			INEEG1.epoch(end+1:end+INEEG2.trials) = INEEG2.epoch(:);
		catch
			disp('pop_mergetset: epoch info removed (information not consistent across datasets)');
		end;
	else
		INEEG1.epoch =[];
	end;

  	% check consistency of merged dataset, regenerate epoch structure
  	% ---------------------------------------------------------------
	if ~isempty(INEEG2.event)
        INEEG1.pnts = size(INEEG1.data,2);
        disp('Reconstituting epoch information...');
        INEEG1 = eeg_checkset(INEEG1, 'eventconsistency');
    end
end

% build the command
% -----------------
if exist('indices') == 1
	com = sprintf('EEG = pop_mergeset( %s, [%s], %d);', inputname(1), int2str(indices), keepall);
else
	com = sprintf('EEG = pop_mergeset( %s, %s, %d);', inputname(1), inputname(2), keepall);		
end

return
