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
% Revision 1.38  2004/02/05 01:37:51  arno
% same
%
% Revision 1.37  2004/02/05 01:29:59  arno
% debug sobi call
%
% Revision 1.36  2004/02/05 01:24:09  arno
% updating sobi calls
%
% Revision 1.35  2003/12/19 02:49:43  arno
% debug rank lowering
%
% Revision 1.34  2003/12/11 16:29:27  arno
% debug history
%
% Revision 1.33  2003/12/11 16:28:01  arno
% remove debug mesg
%
% Revision 1.32  2003/11/07 02:22:18  arno
% fixing history problem
%
% Revision 1.31  2003/10/22 18:10:36  arno
% typo in gui
%
% Revision 1.30  2003/10/15 00:36:33  arno
% same
%
% Revision 1.29  2003/10/15 00:35:41  arno
% debug fastica
%
% Revision 1.28  2003/09/02 21:40:59  lee
% debug last
%
% Revision 1.27  2003/09/02 21:33:52  lee
% changed data rank warning
%
% Revision 1.26  2003/07/26 00:53:32  arno
% computing matrix rank and reducing number of comps accordingly
%
% Revision 1.25  2003/07/24 23:27:24  arno
% debuging ng_ol
%
% Revision 1.24  2003/07/24 00:21:29  arno
% buttons ...
%
% Revision 1.23  2003/07/23 18:36:54  arno
% adding more algos
%
% Revision 1.22  2003/07/22 16:19:32  arno
% tfbss
%
% Revision 1.21  2003/07/22 01:12:22  arno
% adding algos
%
% Revision 1.20  2003/06/29 02:02:30  arno
% last revision backup
%
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

% find available algorithms
% -------------------------
allalgs   = { 'runica' 'binica' 'jader' 'jadeop' 'jade_td_p' 'MatlabshibbsR' 'fastica' ...
              'tica' 'erica' 'simbec' 'unica' 'amuse' 'fobi' 'evd' 'evd24' 'sons' 'sobi' 'ng_ol' ...
              'acsobiro' 'acrsobibpf' 'pearson_ica' 'egld_ica' 'eeA' 'tfbss' 'icaML' 'icaMS' }; % do not use egld_ica => too slow
selectalg = {};
stralg{1} = 'ICA algorithm to use [ ';

linenb    = 1;
count     = 1;
for index = 1:length(allalgs)
    if exist(allalgs{index}) == 2 | exist(allalgs{index}) == 6
        selectalg          = { selectalg{:} allalgs{index} };
        if mod(count+2, 8) == 0, linenb = linenb+1; stralg{linenb} = ''; end;
        if count ~= 1
             stralg{linenb} = [ stralg{linenb} ' | '  allalgs{index} ];
        else stralg{linenb} = [ stralg{linenb}        allalgs{index} ];
        end;
        count               = count+1;
    end;
end;
stralg{linenb} = [ stralg{linenb} ' ] ' ];
    
fig = [];
if nargin < 2 
    % popup window parameters
    % -----------------------
    promptstr    = { strvcat(stralg{:}) ...
                     [ 'Commandline options (See algorithm help messages)' ] };
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
tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean 
switch lower(icatype)
    case 'runica' 
        if nargin < 2
            fig = figure('visible', 'off');
            supergui( fig, {1 1}, [], {'style' 'text' 'string' 'Press Button to interupt runica' }, ...
                      {'style' 'pushbutton' 'string' 'Interupt' 'callback' 'figure(gcbf); set(gcbf, ''tag'', ''stop'');' } );
            drawnow;
        end;
        tmprank = rank(tmpdata(:,1:floor(size(tmpdata,2)/2)));
        if rank(tmpdata) < size(EEG.data,1), 
            disp(['Data rank (' int2str(tmprank) ') less than the number of channels (' int2str(size(EEG.data,1)) ').']);
        end;
        if length(options) < 2
            if rank(tmpdata) == size(EEG.data,1), 
                [EEG.icaweights,EEG.icasphere] = runica( tmpdata, 'lrate', 0.001 );
            else 
                [EEG.icaweights,EEG.icasphere] = runica( tmpdata, 'lrate', 0.001, 'pca', tmprank );
            end;
        else % if there are defined 'options'   
            if rank(tmpdata) == size(EEG.data,1) | ~isempty(findstr('pca', options))
                eval(sprintf('[EEG.icaweights,EEG.icasphere] = runica( tmpdata %s );', options));
            else
                eval(sprintf('[EEG.icaweights,EEG.icasphere] = runica( tmpdata %s, ''pca'', %d );', options, tmprank));
            end;
        end;
     case 'binica'
        if ~isunix | strcmp(computer, 'MAC')
            error('Pop_runica: binica can now only be used under specific UNIX OS');
        end;
        icadefs;
        fprintf(['Warning: if the binary ICA function does not work, check that you have added the\n' ...
                 'binary file location (in the EEGLAB directory) to your Unix /bin directory (.cshrc file)\n']);
        if exist(ICABINARY) ~= 2
            error('Pop_runica: binary ICA program cannot be found. Edit icadefs.m file to specify ICABINARY variable');
        end;
        tmprank = rank(tmpdata(:,1:floor(size(tmpdata,2)/2)));
        if rank(tmpdata) < size(EEG.data,1), 
            disp(['Data rank (' int2str(tmprank) ') less than the number of channels (' int2str(size(EEG.data,1)) ').']);
        end;
        if length(options) < 2
            if rank(tmpdata) == size(EEG.data,1), 
                [EEG.icaweights,EEG.icasphere] = binica( tmpdata, 'lrate', 0.001 );
            else 
                [EEG.icaweights,EEG.icasphere] = binica( tmpdata, 'lrate', 0.001, 'pca', tmprank );
            end;
        else % if defined 'options'
            if rank(tmpdata) == size(EEG.data,1), 
                eval(sprintf('[EEG.icaweights,EEG.icasphere] = binica( tmpdata %s );', options));
            else
                eval(sprintf('[EEG.icaweights,EEG.icasphere] = binica( tmpdata %s, ''pca'', %d );', options, tmprank));
            end;
        end;
     case 'jader' 
        if length(options) < 2
            [EEG.icaweights] = jader( tmpdata );
        else    
            eval(sprintf('[EEG.icaweights] = jader( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'matlabshibbsr' 
        if length(options) < 2
            [EEG.icaweights] = MatlabshibbsR( tmpdata );
        else    
            eval(sprintf('[EEG.icaweights] = MatlabshibbsR( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'pearson_ica' 
        if length(options) < 2
            [tmp EEG.icaweights] = pearson_ica( tmpdata, 'maxNumIterations', 1000,'epsilon',0.0005);
        else    
            eval(sprintf('[tmp EEG.icaweights] = pearson_ica( tmpdata %s );', options));
        end;
        clear tmp;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'egld_ica' 
        disp('Warning: this algorithm is very slow');
        if length(options) < 2
            [tmp EEG.icaweights] = egld_ica( tmpdata );
        else    
            eval(sprintf('[tmp EEG.icaweights] = egld_ica( tmpdata %s );', options));
        end;
        clear tmp;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'eea' 
        if length(options) < 2
            [EEG.icaweights] = eeA( tmpdata );
        else    
            eval(sprintf('[EEG.icaweights] = eeA( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'icaml' 
        if length(options) < 2
            [tmp EEG.icaweights] = icaML( tmpdata );
        else    
            eval(sprintf('[tmp EEG.icaweights] = icaML( tmpdata %s );', options));
        end;
        clear tmp;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'icams' 
        if length(options) < 2
            [tmp EEG.icaweights] = icaMS( tmpdata );
        else    
            eval(sprintf('[tmp EEG.icaweights] = icaMS( tmpdata %s );', options));
        end;
        clear tmp;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'tfbss' 
        if length(options) < 2
             [tmp EEG.icaweights] = tfbss( tmpdata, size(tmpdata,1), 8, 512 );
        else    
            eval(sprintf('[tmp EEG.icaweights] = tfbss( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        clear tmp;
     case 'tica' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icaweights = tica( tmpdata, size(tmpdata,1) );
        else    
            eval(sprintf('EEG.icaweights = tica( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'erica' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icaweights = erica( tmpdata, size(tmpdata,1) );
        else    
            eval(sprintf('EEG.icaweights = erica( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'simbec' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icaweights = simbec( tmpdata, size(tmpdata,1) );
        else    
            eval(sprintf('EEG.icaweights = simbec( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'unica' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icaweights = unica( tmpdata, size(tmpdata,1) );
        else    
            eval(sprintf('EEG.icaweights = unica( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'amuse' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icaweights = amuse( tmpdata );
        else    
            eval(sprintf('EEG.icaweights = amuse( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'fobi' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             [tmp EEG.icaweights] = fobi( tmpdata, size(tmpdata,1) );
        else    
            eval(sprintf('EEG.icaweights = fobi( tmpdata %s );', options));
        end;
        clear tmp;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'evd' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icaweights = evd( tmpdata );
        else    
            eval(sprintf('EEG.icaweights = evd( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'sons' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icaweights = sons( tmpdata );
        else    
            eval(sprintf('EEG.icaweights = sons( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'jadeop' 
        if length(options) < 2
             EEG.icaweights = jadeop( tmpdata );
        else    
            eval(sprintf('EEG.icaweights = jadeop( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'jade_td_p' 
        if length(options) < 2
             EEG.icaweights = jade_td_p( tmpdata, size(tmpdata,1) );
        else    
            eval(sprintf('EEG.icaweights = jade_td_p( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
     case 'evd24' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icaweights = evd24( tmpdata );
        else    
            eval(sprintf('EEG.icaweights = evd24( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'sobi' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        disp('Running sobi...');
        if length(options) < 2
             EEG.icawinv = sobi( EEG.data );
        else    
            eval(sprintf('EEG.icawinv = sobi( EEG.data %s );', options));
        end;
        EEG.icaweights = pinv(EEG.icawinv);
        EEG.icasphere = eye(size(EEG.icaweights,2));
        close(fig);
     case 'ng_ol' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             [tmp EEG.icaweights] = ng_ol( tmpdata );
        else    
            eval(sprintf('[ EEG.icaweights tmp ] = ng_ol( tmpdata %s );', options));
        end;
        EEG.icasphere = eye(size(EEG.icaweights,2));
        clear tmp;
        close(fig);
     case 'acsobiro' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icawinv = acsobiro( EEG.data, size(tmpdata,1) );
        else    
            eval(sprintf('EEG.icawinv = acsobiro( EEG.data %s );', options));
        end;
        EEG.icaweights = pinv(EEG.icawinv);
        EEG.icasphere  = eye(size(EEG.icaweights,2));
        close(fig);
     case 'acrsobibpf' 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        if length(options) < 2
             EEG.icawinv = acrsobibpf( tmpdata, size(tmpdata,1) );
        else    
            eval(sprintf('EEG.icawinv = acrsobibpf( tmpdata %s );', options));
        end;
        EEG.icaweights = pinv(EEG.icawinv);
        EEG.icasphere  = eye(size(EEG.icaweights,2));
        close(fig);
     case 'fastica'
        if ~exist('fastica', 'file') & ~exist('fastica', 'dir')
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
    if length(options) < 2
        com = sprintf('%s = pop_runica(%s, ''%s'');', inputname(1), inputname(1), icatype);
    else
        com = sprintf('%s = pop_runica(%s, ''%s'' %s);', inputname(1),inputname(1), icatype, options);
    end;
end;

return;
