% std_readspec() - returns the stored mean power spectrum for an ICA component 
%                  in a specified dataset.  The spectrum is assumed to have been 
%                  saved in a Matlab file, "[dataset_name].icaspec", in the same
%                  directory as the dataset file. If this file doesn't exist,
%                  use std_spec() to create it or a pre-clustering function
%                  (pop_preclust() or std_preclust()) that calls it. 
% Usage:    
%         >> [spec, freqs] = std_readspec(ALLEEG, setindx, component, freqrange);  
% Inputs:
%   ALLEEG     - a vector of dataset EEG structures (may also be one dataset). 
%                Must contain the dataset of interest (the 'setindx' below).
%   setindx    - [integer] an index of an EEG dataset in the ALLEEG
%                structure for which to read a component spectrum.
%   component  - [integer] index of the component in the selected EEG dataset 
%                for which to return the spectrum
%   freqrange  - [min max in Hz] frequency range to return
%
% Outputs:
%   spec      - the log-power spectrum of the requested ICA component in the
%               specified dataset (in dB)
%   freqs     - vector of spectral frequencies (in Hz)
%
%  See also  std_spec(), pop_preclust(), std_preclust()
%
% Authors:  Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, February, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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
% Revision 1.14  2006/03/28 15:38:13  scott
% help msg
%
% Revision 1.13  2006/03/14 02:32:32  scott
% help msg
%
% Revision 1.12  2006/03/11 07:30:01  arno
% freqrange input
%
% Revision 1.11  2006/03/11 07:25:37  arno
% header
%
% Revision 1.10  2006/03/10 16:33:37  arno
% selecting frequency range for reading
%
% Revision 1.9  2006/03/10 00:37:45  arno
% error msg
%
% Revision 1.8  2006/03/09 18:10:38  arno
% *** empty log message ***
%
% Revision 1.7  2006/03/09 18:10:18  arno
% do not use etc field any more
%
% Revision 1.6  2006/03/09 00:42:09  arno
% fix reading file
%
% Revision 1.5  2006/03/09 00:37:31  arno
% now writing matlab fileend
%
% Revision 1.4  2006/03/09 00:03:57  arno
% read spectrum form matlab file
%
% Revision 1.3  2006/03/08 21:06:37  arno
% rename func
%
% Revision 1.2  2006/03/07 22:21:12  arno
% use fullfile
%

function [X, f] = std_readspec(ALLEEG, abset, comp, freqrange);

X = [];
filename = fullfile( ALLEEG(abset).filepath,[ ALLEEG(abset).filename(1:end-3) 'icaspec' ]);

for k=1:length(comp)
    try,
        specstruct = load( '-mat', filename, [ 'comp' int2str(comp(k)) ], 'freqs' );
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;

    tmpdat    = getfield(specstruct, [ 'comp' int2str(comp(k)) ]);
    if k == 1
        X = zeros(length(comp), length(tmpdat));
    end;
    X(k,:)  = tmpdat;
    f       = specstruct.freqs;
end;

% select frequency range of interest
% ----------------------------------

%check if freqrange is specified
if(~isempty(freqrange))
    maxind = max(find(f <= freqrange(end)));
    minind = min(find(f >= freqrange(1)));
else
    %if not, use whole spectrum
    maxind = length(f);
    minind = 1;
end

f = f(minind:maxind);
X = X(:,minind:maxind);
%if f(end) < freqrange(end)
%    disp(['Warning! Requested high frequency limit, ' ...
%          num2str(freqrange(end)) 'Hz, is out of bounds.']);
%end
%if f(1) > freqrange(1)
%    disp(['Warning! Requested low frequency limit, ' ...
%          num2str(freqrange(1)) 'Hz, is out of bounds.']);
%end

return;
