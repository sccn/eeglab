% eeg_rejmacro() - Internal EEGLAB macro for all pop_ functions that
%                  perform data rejection.   
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
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-08-02 include past rejections in eegplot -ad

% this macro is used by pop_eegplot, pop_eegthresh, pop_const, pop_jointprob, pop_eegthreshfreq, pop_kurtosis

% inputs:
%   rej - input rejection in this workspace (string)
%   rejE - electrode rejection in this workspace (string)
%   elecrange - range of selected electrodes
%   reject - superpose with previous rejections
% optional : nbpnts = number of points per trial

% outputs:
%   rejection arrays ('oldrej', 'oldrejE'), 'command' for eeglot

if ~exist('elecrange')
	help eeg_rejmacro;
	error('Error using eeg_rejmacro');
end;	

if ~exist('nbpnts') nbpnts = EEG.pnts; end;
com2 = [ 'if ~isempty(TMPREJ), [tmprej tmprejE] = eegplot2trial(TMPREJ,' ...
		        int2str(nbpnts) ', EEG.trials, [], []);' ... % include past rejection (see eeg_multieegplot for the color definition
         'tmprejE2 = zeros(EEG.nbchan, length(tmprej));' ...
         sprintf('tmprejE2([%s],:) = tmprejE;', int2str(elecrange)) ... 
         macrorej '= tmprej;' macrorejE '= tmprejE2;' ];
if reject
	if ~exist('icacomp'), icacomp = typerej; end;
    com2 = [com2 sprintf(['%s = pop_rejepoch(%s, tmprej, 1);' ...
		   '[ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET); h(LASTCOM);' ...
	       'eeglab(''redraw''); end;'], inputname(1), inputname(1)); ] ;
else
	com2 = [com2 '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eeglab(''redraw''); end;' ];
end; 

if ~exist('topcommand')
	topcommand = [];
end;

% the first part is used to convert the eegplot output
command = [  com2 topcommand 'clear tmprej tmprejE tmprejE2 TMPREJ;' ];

switch superpose
	case 0, oldrej  =  []; oldrejE =  [];
	case 1, oldrej  = eval(macrorej);
    		oldrejE = eval(macrorejE);
    case 2, oldrej  = EEG.reject.rejglobal;
    		oldrejE = EEG.reject.rejglobalE;
end;

if ~isempty(oldrejE) oldrejE = oldrejE(elecrange,:); end;

