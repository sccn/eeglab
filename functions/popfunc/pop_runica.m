% pop_runica() - Run an ICA decomposition on an EEG dataset 
%                using runica(),binica(), or other ICA algorithm.
% Usage:
%   >> OUT_EEG = pop_runica( IN_EEG ); % pops-up a data entry window
%   >> OUT_EEG = pop_runica( IN_EEG, 'key', 'val' ); % no pop_up
%
% Graphic interface:
%   "ICA algorithm to use" - [edit box] The type of ICA algorithm 
%                 to use for the ICA decomposition. Command line
%                 equivalent: 'icatype'
%   "Commandline options" - [edit box] Command line options to forward
%                 to the ICA algorithm. Command line equivalent: 'options' 
% Inputs:
%   IN_EEG      - input EEG dataset
%
% Optional inputs:
%   'icatype'   - ['runica'|'binica'|'jader'|'fastica'] ICA algorithm 
%                 to use for the ICA decomposition. The nature of any 
%                 differences in the results of these algorithms have 
%                 not been well characterized. Default is binica(), if
%                 found, else runica().
%   'dataset'   - [integer array] dataset indices.
%   'key','val' - ICA algorithm options (see ICA routine help messages).
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
% Revision 1.54  2005/07/26 23:55:17  arno
% rewrote the function; allow to process several datasets
%
% Revision 1.53  2005/03/14 22:20:45  arno
% fixing hilit correction
%
% Revision 1.52  2005/03/14 22:13:59  hilit
% fixing a typo
%
% Revision 1.51  2005/03/14 19:43:28  arno
% saving old weights
%
% Revision 1.50  2005/03/13 19:38:10  scott
% saving EEG.oldwts and EEG.oldsph as cell arrays containing all former wts and sph
%
% Revision 1.49  2005/03/13 17:45:04  peter
% saving wts and sph as EEG.oldwts (cell array) and EEG.oldsph
% before clearing them
% NOTE: 'binica' with 'maxsteps',5  doesnt work!!! (maxsteps not changed) ??
%
% Revision 1.48  2005/03/07 19:53:42  arno
% fixing fastica call
%
% Revision 1.47  2005/02/02 01:29:28  arno
% new method to compute weights, icawinv...
%
% Revision 1.46  2005/02/02 01:21:40  arno
% remove old ICA information
%
% Revision 1.45  2005/02/01 23:09:42  arno
% call to icaML and icaMS
%
% Revision 1.44  2004/09/09 00:13:31  arno
% fixing rank problem
%
% Revision 1.43  2004/08/20 18:43:07  arno
% remove MAC error for Osx
%
% Revision 1.42  2004/08/20 18:11:11  arno
% no eval for binica
%
% Revision 1.41  2004/05/20 15:48:19  arno
% debug fastica call
%
% Revision 1.40  2004/03/22 22:15:39  arno
% pca option for binica
%
% Revision 1.39  2004/02/05 01:50:52  arno
% call for acsobiro
%
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

function [ALLEEG, com] = pop_runica( ALLEEG, varargin )

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

linenb    = 1;
count     = 1;
for index = length(allalgs):-1:1
    if exist(allalgs{index}) ~= 2 && exist(allalgs{index}) ~= 6
        allalgs(index) = [];
    end;
end;

% decode input arguments
% ----------------------
pop_up = 0;
if nargin < 2 
    pop_up = 1;
    g.defaultdataset = 1;
else
    if ~isstr(varargin{1})
        g.icatype = varargin{1};
        options    = varargin(2:end);
    else
        [ g options ] = finputcheck( varargin, { 'icatype'        'string'  allalgs   'runica'; ...
                                                 'dataset'        'integer' []        1;
                                                 'defaultdataset' 'integer' []        [] }, ...
                                 'pop_runica', 'ignore');
        if isstr(g), error(g); end;
        if ~isempty(g.defaultdataset), pop_up = 1; end;
    end;    
end;

fig = [];
if pop_up
    % popup window parameters
    % -----------------------
    promptstr    = { { 'style' 'text'    'string' 'ICA algorithm to use (click to select)' } ...
                     { 'style' 'listbox' 'string' strvcat(allalgs{:}) } ...
                     { 'style' 'text'    'string' 'Commandline options (See algorithm help messages)' } ...
                     { 'style' 'edit'    'string' '' }  };
	inistr       = { 'runica' '' };
    result       = inputgui( { [2 1] [2 1] }, promptstr, 'pophelp(''pop_runica'')', ...
                             'Run ICA decomposition -- pop_runica()');
    if length(result) == 0 return; end;
    g.icatype      = allalgs{result{1}};
    g.dataset      = 1:length(ALLEEG);
    options        = eval( [ '{' result{2} '}' ]);
end;

% select datasets, create new big dataset if necessary
% ----------------------------------------------------
if length(g.dataset) == 1
    EEG = ALLEEG(g.dataset)
else
    disp('Concatenating datasets...');
    EEG = ALLEEG(g.dataset(1));
    
    % compute total data size
    % -----------------------
    totalpnts = 0;
    for i = g.dataset
        totalpnts = totalpnts+ALLEEG(g.dataset(i)).pnts*ALLEEG(g.dataset(i)).trials;
    end;
    EEG.data = zeros(EEG.nbchan, totalpnts);
    
    % copy data
    % ---------
    cpnts = 1;
    for i = g.dataset
        tmplen = ALLEEG(g.dataset(i)).pnts*ALLEEG(g.dataset(i)).trials;
        EEG.data(:,cpnts:cpnts+tmplen-1) = ALLEEG(g.dataset(i)).data(:,:);
        cpnts = cpnts+1;
    end;
    EEG.icaweights = [];
    EEG.trials = 1;
    EEG.pnts   = size(EEG.data,2);
end;    

% Store and then remove current EEG ICA weights and sphere
% ---------------------------------------------------
fprintf('\n');
if ~isempty(EEG.icaweights)
    fprintf('Saving current ICA weights in "EEG.etc" sub-structure.\n');
    if ~isfield(EEG,'etc'), EEG.etc = []; end;
    if ~isfield(EEG.etc,'oldicaweights')
        EEG.etc.oldicaweights = {};
        EEG.etc.oldicasphere = {};
    end;
    EEG.etc.oldicaweights = { EEG.icaweights EEG.etc.oldicaweights{:} };
    EEG.etc.oldicasphere  = { EEG.icasphere  EEG.etc.oldicasphere{:}  };
end
EEG.icaweights = [];
EEG.icasphere  = [];
EEG.icawinv    = [];
EEG.icaact     = [];

% is pca already an option?
% -------------------------
pca_opt = 0;
for i = 1:length(options)
    if isstr(options{i})
        if strcmpi(options{1}, 'pca')
            pca_opt = 1;
        end;
    end;
end;

%------------------------------
% compute ICA on a definite set
% -----------------------------
tmpdata = reshape( EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean 
switch lower(g.icatype)
    case 'runica' 
        if nargin < 2
            fig = figure('visible', 'off');
            supergui( fig, {1 1}, [], {'style' 'text' 'string' 'Press Button to interupt runica' }, ...
                      {'style' 'pushbutton' 'string' 'Interupt' 'callback' 'figure(gcbf); set(gcbf, ''tag'', ''stop'');' } );
            drawnow;
        end;
        tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2))));
        if tmprank == size(EEG.data,1) | pca_opt
            [EEG.icaweights,EEG.icasphere] = runica( tmpdata, 'lrate', 0.001, options{:} );
        else 
            disp(['Data rank (' int2str(tmprank) ') less than the number of channels (' int2str(size(EEG.data,1)) ').']);
            [EEG.icaweights,EEG.icasphere] = runica( tmpdata, 'lrate', 0.001, 'pca', tmprank, options{:} );
        end;
     case 'binica'
        icadefs;
        fprintf(['Warning: IF the binary ICA function does not work, check that you have added the\n' ...
                 'binary file location (in the EEGLAB directory) to your Unix /bin directory (.cshrc file)\n']);
        if exist(ICABINARY) ~= 2
            error('Pop_runica: binary ICA program cannot be found. Edit icadefs.m file to specify ICABINARY variable');
        end;
        tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2))));
        if tmprank == size(EEG.data,1) | pca_opt
            [EEG.icaweights,EEG.icasphere] = binica( tmpdata, 'lrate', 0.001, options{:} );
        else 
            disp(['Data rank (' int2str(tmprank) ') less than the number of channels (' int2str(size(EEG.data,1)) ').']);
            [EEG.icaweights,EEG.icasphere] = binica( tmpdata, 'lrate', 0.001, 'pca', tmprank, options{:} );
        end;
     case 'pearson_ica' 
        if isempty(options)
            disp('Warning: EEG default for pearson ICA changed to 1000 iterations and epsilon=0.0005');
            [tmp EEG.icaweights] = pearson_ica( tmpdata, 'maxNumIterations', 1000,'epsilon',0.0005);
        else    
            [tmp EEG.icaweights] = pearson_ica( tmpdata, options{:});
        end;
     case 'egld_ica', disp('Warning: this algorithm is very slow!!!');
                      [tmp EEG.icaweights] = egld_ica( tmpdata, options{:} );
     case 'tfbss' 
        if  isempty(options)
             [tmp EEG.icaweights] = tfbss( tmpdata, size(tmpdata,1), 8, 512 );
        else    
             [tmp EEG.icaweights] = tfbss( tmpdata, options{:} );
        end;
     case 'jader',         [EEG.icaweights] = jader( tmpdata, options{:} );
     case 'matlabshibbsr', [EEG.icaweights] = MatlabshibbsR( tmpdata, options{:} );
     case 'eea',           [EEG.icaweights] = eeA( tmpdata, options{:} );
     case 'icaml',         [tmp EEG.icawinv] = icaML( tmpdata, options{:} );
     case 'icams',         [tmp EEG.icawinv] = icaMS( tmpdata, options{:} );
     case 'fastica',       [ ICAcomp, EEG.icawinv, EEG.icaweights] = fastica( tmpdata, 'displayMode', 'off', options{:} );
     case { 'tica' 'erica' 'simbec' 'unica' 'amuse' 'fobi' 'evd' 'sons' ...
            'jadeop' 'jade_td_p' 'evd24' 'sobi' 'ng_ol' 'acrsobiro' 'acrsobibpf' } 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        
        if isempty(options), options = { size(tmpdata,1) }; end;
        switch lower(g.icatype)
         case 'tica',     EEG.icaweights = tica( tmpdata, options{:} );
         case 'erica',    EEG.icaweights = erica( tmpdata, options{:} );
         case 'simbec',   EEG.icaweights = simbec( tmpdata, options{:} );
         case 'unica',    EEG.icaweights = unica( tmpdata, options{:} );
         case 'amuse',    EEG.icaweights = amuse( tmpdata );
         case 'fobi',     [tmp EEG.icaweights] = fobi( tmpdata, options{:} );
         case 'evd',      EEG.icaweights = evd( tmpdata, options{:} );
         case 'sons',     EEG.icaweights = sons( tmpdata, options{:} );
         case 'jadeop',   EEG.icaweights = jadeop( tmpdata, options{:} );
         case 'jade_td_p',EEG.icaweights = jade_td_p( tmpdata, options{:} );
         case 'evd24',    EEG.icaweights = evd24( tmpdata, options{:} );
         case 'sobi',     EEG.icawinv = sobi( EEG.data, options{:} );
         case 'ng_ol',    [tmp EEG.icaweights] = ng_ol( tmpdata, options{:} );
         case 'acsobiro', EEG.icawinv = acsobiro( EEG.data, options{:} );
         case 'acrsobibpf', EEG.icawinv = acrsobibpf( tmpdata, options{:} );
        end;
        clear tmp;
        close(fig);
     otherwise, error('Pop_runica: unrecognized algorithm');
end;

% update weight and inverse matrices etc...
% -----------------------------------------
if ~isempty(fig), try, close(fig); catch, end; end;
if isempty(EEG.icasphere)
    EEG.icasphere  = eye(size(EEG.icaweights,2));
end;
if isempty(EEG.icaweights)
    EEG.icaweights = pinv(EEG.icawinv);
end;
if isempty(EEG.icawinv)
    EEG.icawinv    = pinv(EEG.icaweights*EEG.icasphere); % a priori same result as inv
end;

% copy back data to datasets if necessary
% ---------------------------------------
if length(g.dataset) > 1
    for i = 1:g.dataset
        ALLEEG(i).icaweights = EEG.icaweights;
        ALLEEG(i).icasphere  = EEG.icasphere;
        ALLEEG(i).icawinv    = EEG.icawinv;
    end;            
    ALLEEG = eeg_checkset(ALLEEG);
else
    EEG = eeg_checkset(EEG);
    ALLEEG(g.dataset) = EEG;
end;

if pop_up
    com = sprintf('%s = pop_runica(%s, %s);', inputname(1),inputname(1), ...
                  vararg2str({ 'icatype' g.icatype 'dataset' g.dataset options{:} }) );
end;

return;
