% eeg_retrieve() - Retrieve an EEG dataset from the variable
%                  containing all datasets (standard: ALLEEG).
%
% Usage: >> EEG = eeg_retrieve( ALLEEG, index );
%
% Inputs:
%   ALLEEG     - variable containing all datasets
%   index      - index of the dataset to retrieve
%
% Outputs:
%   EEG        - output dataset. The workspace variable EEG is also updated 
%   ALLEEG     - updated ALLEEG structure
%   CURRENTSET - workspace variable index of the current dataset
%
% Note: The function performs >> EEG = ALLEEG(index);
%       It also runs eeg_checkset() on it.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeg_store(), eeglab()

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
% Revision 1.15  2005/09/27 22:06:53  arno
% fix empty dataset retrieve
%
% Revision 1.14  2005/09/08 21:58:58  arno
% preserving 'saved' field
%
% Revision 1.13  2005/09/08 21:37:31  arno
% fix retrieving multiple datasets
%
% Revision 1.12  2005/09/08 16:54:29  arno
% preserve saved field
%
% Revision 1.11  2005/09/08 16:41:04  arno
% fixing storing problem for saved datasets
%
% Revision 1.10  2005/08/08 18:03:55  arno
% fix dissimilar structure problem
%
% Revision 1.9  2005/08/04 23:37:24  arno
% typo
%
% Revision 1.8  2005/08/02 16:54:15  arno
% adding outputs
%
% Revision 1.7  2005/08/01 22:43:26  arno
% eeg_options call
%
% Revision 1.6  2005/08/01 15:44:31  arno
% do not load data when retrieveing several datasets
%
% Revision 1.5  2005/08/01 14:50:44  arno
% loaddata
%
% Revision 1.4  2002/08/11 17:32:00  arno
% header
%
% Revision 1.3  2002/07/22 21:06:52  arno
% nothing
%
% Revision 1.2  2002/04/23 23:48:21  arno
% standalone version
%
% Revision 1.1  2002/04/18 19:49:29  arno
% Initial revision
%
% Revision 1.1  2002/04/05 17:46:04  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 

function [EEG, ALLEEG, CURRENTSET] = eeg_retrieve( ALLEEG, CURRENTSET);

if nargin < 2
	help eeg_retrieve;
	return;
end;	

%try
    eeg_optionsbackup;
    eeg_options;

    try,   % check whether recent changes to this dataset have been saved or not
           %--------------------------------------------------------------------
           tmpsaved = { ALLEEG(CURRENTSET).saved };
    catch, tmpsaved = 'no';
    end;

    if length(CURRENTSET) > 1 & option_storedisk
        [ EEG tmpcom ] = eeg_checkset(ALLEEG(CURRENTSET)); % do not load data if several datasets
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    else
        if CURRENTSET ~= 0
            [ EEG tmpcom ] = eeg_checkset(ALLEEG(CURRENTSET), 'loaddata');
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        else
            EEG = eeg_emptyset; % empty dataset
            return;
        end;
    end;
    
    % retain saved field
    % ------------------
    for index = 1:length(CURRENTSET)
        ALLEEG(CURRENTSET(index)).saved = tmpsaved{index};
        EEG(index).saved                = tmpsaved{index};
    end;
    
%catch
%	fprintf('Warning: cannot retrieve dataset with index %d\n', CURRENTSET); 
%	return;
%end;

return;

