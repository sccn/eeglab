% pop_runica() - Run an ICA decomposition on an EEG dataset 
%                using runica(),binica(), or other ICA algorithm.
% Usage:
%   >> OUT_EEG = pop_runica( IN_EEG ); % pops-up a data entry window
%   >> OUT_EEG = pop_runica( IN_EEG, ica_type, options ); % no pop_up
%
% Graphic interface:
%   "ICA algorithm to use" - [edit box] The type of ICA algorithm 
%                 to use for the ICA decomposition. 
%                 equivalent: 'rhe ica_type'
%   "Commandline options" - [edit box] Command line options to forward
%                 to the ICA algorithm. Command line eqivalent: 'options' 
% Inputs:
%   IN_EEG      - input EEG dataset
%   ica_type    - ['runica'|'binica'|'jader'|'fastica'] ICA algorithm 
%                 to use for the ICA decomposition. The nature of any 
%                 differences in the results of these algorithms have 
%                 not been well characterized. Default is binica(), if
%                 found, else runica().
%   options     - ICA algorithm options (see ICA routine help messages).
% 
% Note:
% 1) Infomax is the ICA algorithm we use most. It is based on Tony Bell's
%    algorithm implemented for automated use by Scott Makeig using the 
%    natural gradient of Amari et al.. It can also extract sub-Gaussian 
%    sources using the 'extended' ICA option of Lee and Girolami. Function
%    runica() is the all-Matlab version; binica() calls the (1.5x faster) 
%    binary version (separate download) translated to C by Sigurd Enghoff.
% 2) jader() calls the JADE algorithm of Jean-Francois Cardoso
%    It is included in the EEGLAB toolbox by his permission. 
%    See >> help jader
% 3) To run fastica(), download the fastICA toolbox from
%    http://www.cis.hut.fi/projects/ica/fastica/ and make it available 
%    in your Matlab path. According to the authors, default parameters
%    are not optimal: Try 'approach', 'sym' to estimate components in
%    parallel.
%
% Outputs:
%   OUT_EEG = Input EEGLAB dataset with new .weights and .sphere field values.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: runica(), binica(), jader(), fastica()

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
% Revision 1.18  2003/02/23 08:38:32  scott
% header edit -sm
%
% Revision 1.17  2003/02/19 19:22:40  arno
% updating header for GUI
%
% Revision 1.16  2003/01/15 22:07:59  arno
% typo
%
% Revision 1.15  2002/12/24 01:27:55  arno
% debug for 'pca' option
%
% Revision 1.14  2002/12/05 03:12:06  arno
% fixing fig problem
%
% Revision 1.13  2002/11/15 18:01:06  arno
% adding more warning messages
%
% Revision 1.12  2002/11/13 23:05:53  arno
% problem from command line call
%
% Revision 1.11  2002/11/13 17:08:26  scott
% help msg
% .,
%
% Revision 1.10  2002/10/25 23:55:09  arno
% interupt only for runica
%
% Revision 1.9  2002/10/25 23:52:01  arno
% debugging for Mac
%
% Revision 1.8  2002/10/23 18:09:58  arno
% new interupt button
%
% Revision 1.7  2002/08/23 15:04:29  scott
% help msg
%
% Revision 1.6  2002/08/19 21:53:40  arno
% same
%
% Revision 1.5  2002/08/19 21:37:29  arno
% debugging fastica for mac
%
% Revision 1.4  2002/08/12 02:28:24  arno
% inputdlg2
%
% Revision 1.3  2002/05/02 21:39:42  arno
% editing message
%
% Revision 1.2  2002/04/18 16:01:24  scott
% editted msgs -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 add the eeglab options -ad
% 03-18-02 add other decomposition options -ad
% 03-19-02 text edition -sm

function [EEG, com] = pop_runica( EEG, icatype, varargin )

com = '';
if nargin < 1   
    help pop_runica;
    return;
end;

fig = [];
if nargin < 2 
    % popup window parameters
    % -----------------------
    promptstr    = { [ 'ICA algorithm to use [ runica | binica | jader | fastICA ]' ] ...
                      ['Commandline options (See algorithm help messages)']};
	inistr       = { 'runica' '' };
	result       = inputdlg2( promptstr, 'Run ICA decomposition -- pop_runica()', 1,  inistr, 'pop_runica');
	if length(result) == 0 return; end;
	icatype      = result{1};
	options      = [ ',' result{2} ];
else
	options = [];
	for i=1:length( varargin )
		if isstr( varargin{ i } )
			options = [ options ', ''' varargin{i} '''' ];
		else
			options = [ options ', [' num2str(varargin{i}) ']' ];
		end;
	end;	
end;

%------------------------------
% compute ICA on a definite set
% -----------------------------
tmpdata = reshape( EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
switch lower(icatype)
    case 'runica' 
        if nargin < 2
            fig = figure('visible', 'off');
            supergui( fig, {1 1}, [], {'style' 'text' 'string' 'Press Button to interrupt runica' }, ...
                      {'style' 'pushbutton' 'string' 'Interupt' 'callback' 'figure(gcbf); set(gcbf, ''tag'', ''stop'');' } );
            drawnow;
        end;
        if length(options) < 2
            [EEG.icaweights,EEG.icasphere] = runica( tmpdata, 'lrate', 0.001 );
        else    
            eval(sprintf('[EEG.icaweights,EEG.icasphere] = runica( tmpdata %s );', options));
        end;
     case 'binica'
        if ~isunix | strcmp(computer, 'MAC')
            error('Pop_runica: binica can now only be used under specific UNIX OS');
        end;
        icadefs;
        fprintf(['Warning: if the binary function does not work, check that you have added the\n' ...
                 'binary file location (in the eeglab directory) to you BIN Unix directory (.cshrc file)\n']);
        if exist(ICABINARY) ~= 2
            error('Pop_runica: binary ica program cannot be found. Edit icadefs.m file to specify ICABINARY variable');
        end;
        if length(options) < 2
            [EEG.icaweights,EEG.icasphere] = binica( tmpdata, 'lrate', 0.001 );
        else    
            eval(sprintf('[EEG.icaweights,EEG.icasphere] = binica(tmpdata %s );', options));
        end;
     case 'jader' 
        if length(options) < 2
            [EEG.icaweights] = jader( tmpdata );
        else    
            eval(sprintf('[EEG.icaweights] = jader( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'fastica'
        if exist('fastica') ~= 2
            error('Pop_runica: to use fastica, you must first download the toolbox (see >> help pop_runica)');
        end;     
        if length(options) < 2
            eval([ '[ ICAcomp, EEG.icaweights,EEG.icasphere] = fastica( tmpdata, ''displayMode'', ''off'' );' ]);
        else    
            eval(sprintf('[ ICAcomp, EEG.icaweights,EEG.icasphere] = fastica( tmpdata, ''displayMode'', ''off'' %s );', options));
        end;
     otherwise, error('Pop_runica: unrecognized algorithm');
end;
if ~isempty(fig), try, close(fig); catch, end; end;
EEG.icawinv    = pinv(EEG.icaweights*EEG.icasphere); % a priori same result as inv

eeg_options; 
if option_computeica
    EEG.icaact    = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
    EEG.icaact    = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end;
if nargin < 2
    if length(options < 2)
        com = sprintf('%s = pop_runica(%s, ''%s'');', inputname(1), inputname(1), icatype);
    else
        com = sprintf('%s = pop_runica(%s, ''%s'' %s);', inputname(1), icatype, options);
    end;
end;
return;
