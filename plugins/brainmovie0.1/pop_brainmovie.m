% pop_brainmovie() - function to make movie from an EEGLAB dataset
% 
% Usage:
%  pop_brainmovie(ALLEEG, 'key', 'val', ...);
%
% Inputs:
%    ALLEEG    - array of EEG datasets with 2 conditions
%
% Optional inputs:
%  'mode'     - ['compute'|'movie'|'computemovie'|'auto'] the function can take
%               a lot of time to compute the time-frequency decomposition.
%               moreover, you only want to compute it once and then adjst
%               the parameter of your movie(s). 
%               'compute' only performs the time-frequency decompostion and
%                         save the results on disk
%               'movie' only load the time-frequency decompostion on disk and 
%                       generate a movie.
%               'computemovie' does both. 
%               Default is 'auto' if there is no time-frequency 
%               decomposition in the time-freq, use 'computemovie' otherwise
%               uses 'movie' mode.
%  'comps'      - [integer vector] integer vector for computing the time-
%               decompotion of ICA components and to include in the movie.
%               By defaults, all the components are used.
%  'freqparams' - [cell array] time-frequency decomposition parameters for
%               timef and crossf. The first parameter must be the number of
%               cycles. Default is { 3, 'winsize' 96, 'padratio', 4, 'maxfreq', 25,
%               'alpha', 0.01, 'subitc', 'on' }.
%  'diffmovie'  - ['on'|'off'] plot the difference movie or simply the movie for
%               the given datasets. For the difference movies, there must only
%               be two input datasets (ALLEEG length must be equal to 2). The movie
%               is divided in 3 rows, the first row depict the first condition, the
%               second row the second condition and the third row the difference.
%               Default is 'off'.
%  'freqs'      - [real vector] of frequency values in Hz. For each frequency value
%               the closest value in the frequency decomposition will be found and
%               a movie will be generated. Default is all frequencies.
%  'quality'    - ['ultrafast'|'fast'|'getframe'|'slow'] different output speed and
%               quality for generating movies. Default is 'ultrafast' where the function
%               does not call makemovie(). For ['fast'|'getframe'|'slow'], see
%               makemovie function help.
%  'confirm'    - ['on'|'off'] ask for confimartion for time-consuming computation.
%               Default is 'on'.
%
% Name and folder inputs:
%  'moviename'   - [string] movie name. If you used the frequency input
%                to generate movies at several frequencies, the frequency
%                value will be added at the end of this name for each
%                movie. Default is 'movie'.
%  'moviefolder' - [string] output folder for the movie. Default is the current 
%                folder.
%  'tfname'      - [string] used a a basename to load and save infos.
%                Default is 'tfparams'.
%  'tffolder'    - [string] folder name were the time-frequency decomposition need to
%                be saved/read. Default is the current directory.
%  'framefolder' - [string] output folder for the frames (in .fig format). Default 
%                is a folder named 'movieframes' creqated in the current folder.
%
% Brainmovie arguments:
%  'showcomps'   - [integer vector] integer vector for showing independent
%                components. Default is comps.
%  'coordinates' - [real array] component equivalent dipole locations. Size is
%                (nbdipoles,2). If this argument is not defined, components are
%                plotted on a circle. If this argument is not defined and the 'sources'
%                field of EEG is defined, the function uses the BESA dipole locations.
%  'circfactor'  - [real array] brainmovie 'circfactor' argument.
%  'movparams'   - [cell array or string] brainmovie parameters. If a string is given
%                as argument, 'mriside' will plot a brainmovie using an average
%                structural side MRI image as background, 'mritop' will plot a 
%                brainmovie from the top and 'mrirear' will plot plot a brainmovie from 
%                the rear. Default is 'mriside'.
%  'addmovparams' - [cell array] brainmovie additional parameters. Allow to use the
%                standard movie parameters and change some of the parameters.
%  'eventprob'   - [event] include response probability of events of type given as
%                a input. The event type name must be present in the EEG dataset.
%                By defining for instance 'eventprob', 'rt' it is possible to plot the 
%                probability of response of the subject assumning events of type 'rt'
%                are present in the dataset.
%  'oneframe'    - ['on'|'off'] only plot one slide of the movie (the first one). This
%                can be usefull to optimize the movie parameters.
%
% Makemovie arguments:
%  'makemovie'   - [cell array]  makemovie() parameters. See >> help makemovie
%                for makemovie() arguments (movie output quality, crop ...). 
%                Default is output with same size as the plotted figure and using
%                { 'mode' 'fast' }.
%
% Example:
%   pop_brainmovie(ALLEEG, 'comps', [1 3 4 5 6], 'freqs', 5);
%   % compute time-frequency decomposition and then the movie at 6=5 Hz assuming 
%   % the datasets have dipoles equivalent for its components.
%
%   pop_brainmovie(ALLEEG, 'comps', [1 3 4 5 6], 'freqs', 5, 'eventprob', 'rt');
%   % same as above but add the reaction time probability of response
%   % note that the time-frequency decomposition does not have to be recomputed
%   % since the previous time-frequency decomposition is in the current folder.
%   % This can be tuned manually using the 'mode' parameter.
%
% Note: 1) do not forget to specify the 'tfname' and 'moviename' parameters if you
%       intend to generate several movies in the same folder (otherwise the time
%       frequency decomposition and the movie itself will be overwritten).
%       2) always run the pop_brainmovie command in the directory were the files 
%       were saved. Otherwise set the 'tffolder' parameter.
%
% Author: Arnaud Delorme
%
% See also: brainmovie(), timecrossf()

% $Log: not supported by cvs2svn $
% Revision 1.22  2002/12/06 19:08:19  arno
% plotorder debug
%
% Revision 1.21  2002/12/06 18:47:10  arno
% plotting movie from the rear
%
% Revision 1.20  2002/12/05 15:49:42  arno
% debuging plotorder
%
% Revision 1.19  2002/12/04 22:50:22  arno
% implementing plotorder
%
% Revision 1.18  2002/11/26 18:58:20  arno
% scale power duplicate debugging
%
% Revision 1.17  2002/11/26 18:42:50  arno
% moviefolder now functionnal
%
% Revision 1.16  2002/11/26 16:16:51  arno
% always including frequency in movie name
% ,
%
% Revision 1.15  2002/11/22 02:10:10  arno
% still commenting the revertangle function
%
% Revision 1.14  2002/11/21 22:57:17  arno
% commenting revertangle
%
% Revision 1.13  2002/11/21 20:29:38  arno
% positive and neg scale
%
% Revision 1.12  2002/11/21 19:06:22  arno
% more debugging
%
% Revision 1.11  2002/11/21 18:04:38  arno
% deleting previous tf data
%
% Revision 1.10  2002/11/21 02:08:18  arno
% more debugging of freqs ...
%
% Revision 1.9  2002/11/20 19:22:10  arno
% more besa and allenv debugging
%
% Revision 1.8  2002/11/20 19:08:52  arno
% searching for besa dipole locations
%
% Revision 1.7  2002/11/20 15:33:24  arno
% more correction, some thanks to cooper
%
% Revision 1.6  2002/11/20 01:42:12  arno
% debugging
%
% Revision 1.5  2002/11/20 01:10:48  arno
% header update
%
% Revision 1.4  2002/11/20 00:50:58  arno
% testing
%
% Revision 1.3  2002/11/20 00:46:44  arno
% default subitc
%
% Revision 1.2  2002/11/19 23:15:41  arno
% debugging
%
% Revision 1.1  2002/11/18 23:16:04  arno
% Initial revision
%
% Revision 1.2  2002/11/07 21:55:28  arno
% new default for timesout
%
% Revision 1.1  2002/08/02 16:07:11  arno
% Initial revision
%

function pop_brainmovie(ALLEEG, varargin);

if nargin < 2
	help pop_brainmovie;
	return;
end;
g = finputcheck(varargin, { 'mode'	      'string'        { 'compute' 'movie' 'computemovie' 'auto' }     'auto';
                            'comps'       'integer'       [1 Inf]                                  1:size(ALLEEG(1).icaact,1);
							'freqparams'  'cell'          {}                                       {};
							'diffmovie'   'string'        { 'on' 'off' }                           'off';
							'confirm'     'string'        { 'on' 'off' }                           'on';
							'moviename'   'string'		  {}									   'movie';
							'moviefolder' 'string'        {}                                       '';
							'tfname'      'string'		  {}									   'tfparams';
							'tffolder'    'string'        {}                                       '';
							'framefolder' 'string'        {}                                       [ addfinalsep(pwd) 'movieframes'];
							'movparams'   {'string' 'cell'}       []							   'mriside';
                            'addmovparams'  'cell'          {}							           {};
							'showcomps'   'integer'       []									   [];
							'coordinates' 'real'          []                                       [];
                            'circfactor'  'real'          []                                       [];
                            'freqs'       'real'          []                                       [];
                            'oneframe'    'string'        { 'on' 'off' }                           'off';
                            'quality'     'string'        { 'ultrafast' 'fast' 'getframe' 'slow' } 'ultrafast';
                            'makemovie'   'cell'          {}                                       {};
                            'eventprob'   ''              []                                       [] });
if isstr(g), error(g); end;
clear functions;

% checking parameters
% -------------------
if strcmpi(g.diffmovie, 'on') &	length(ALLEEG) ~= 2
    error('For difference movies: nedd exactly to process 2 datasets');
end;

% cfunding components to show
% ---------------------------
if isempty(g.showcomps), g.showcomps = g.comps; end;
for index = 1:length(g.showcomps)
    tmpshow(index) = find(g.showcomps(index) ==  g.comps);
    if isempty(tmpshow(index))
        error('showcomps error: component not found');
    end;
end;
g.showcomps = tmpshow;

if length(unique(cell2mat({ALLEEG(:).pnts}))) > 1
    error('All datasets must have the same number of points');
end;
if isempty(g.freqs) & strcmp(g.confirm, 'on')
    disp('********** USER ATTENTION REQUIRED ************');
    r = input('Are you sure you want to make a movie a each of the output frequencies (y/n):', 's');
    if r(1) == 'n', disp('Cancelling movie call'); return; end;
end;
g.tffolder    = addfinalsep(g.tffolder);
g.moviefolder = addfinalsep(g.moviefolder);
g.framefolder = addfinalsep(g.framefolder);

% spectral options
% ----------------
if isempty(g.freqparams)
    g.freqparams = { 3, 'winsize', 96, ...
                     'timesout' , min(200, ALLEEG(1).pnts), ...
                     'padratio', 4, ...
                     'maxfreq', 25, ...
                     'alpha', 0.01, ...
                     'subitc', 'on' };
end;
g.freqparams = { ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, ALLEEG(1).srate, ...
                 g.freqparams{:}, ...
                 'plotersp', 'off', ...
                 'plotitc', 'off', ...
                 'plotamp', 'off', ...
                 'plotphase', 'off' };

% read all parameters
% -------------------
nbconditions = length(ALLEEG);

if strcmpi(g.mode, 'compute') | strcmpi(g.mode, 'computemovie') | ...
        ( strcmpi(g.mode, 'auto') & ~exist([g.tffolder g.tfname '_freqs']))
    if strcmpi(g.mode, 'auto')
        fprintf('Auto mode: %s file does not exist, running time-freq. decomposition\n', [g.tffolder g.tfname '_freqs']);
    end;
    if strcmp(g.confirm, 'on')
        disp('********** USER ATTENTION REQUIRED ************');
        r = input('Are you sure you want to compute time-frequency decompositions (y/n):', 's');
        if r(1) == 'n', disp('Cancelling movie call'); return; end;
    end;
    
	% compute timef and crossf for all components
	% -------------------------------------------
	if strcmpi(g.diffmovie, 'on')
        [ALLERSP ALLITC ALLCROSSF ALLCROSSFANGLE times freqs] ...
			= timecrossf( { ALLEEG(1).icaact(g.comps,:) ALLEEG(2).icaact(g.comps,:)} , g.freqparams{:});
	else 
        for ind =1:nbconditions
            [ ALLERSP(:,ind), ALLITC(:,ind), ALLCROSSF(:,:,ind), ALLCROSSFANGLE(:,:,ind), times, freqs] ...
                = timecrossf( ALLEEG(ind).icaact(g.comps,:), g.freqparams{:});
        end;
    end;
	
	eval(['save ' g.tffolder g.tfname '_freqs freqs']);
	eval(['save ' g.tffolder g.tfname '_ALLERSP ALLERSP']);
	eval(['save ' g.tffolder g.tfname '_ALLITC ALLITC']);
	eval(['save ' g.tffolder g.tfname '_ALLCROSSF ALLCROSSF']);
	eval(['save ' g.tffolder g.tfname '_ALLCROSSFANGLE ALLCROSSFANGLE']);
	eval(['save ' g.tffolder g.tfname '_times times']);
	eval(['save ' g.tffolder g.tfname '_freqs freqs']);
	disp('**************** Computation terminated and saved');
    system(['rm -f ' g.tffolder g.tfname '_newERSP' ]);
else
    if strcmpi(g.mode, 'auto')
        fprintf('Auto mode: found existing files\n');
    end;
	disp('Loading files');
	eval(['load ' g.tffolder g.tfname '_ALLERSP' ]); 
	eval(['load ' g.tffolder g.tfname '_ALLITC' ]); 
	eval(['load ' g.tffolder g.tfname '_ALLCROSSF' ]); 
	eval(['load ' g.tffolder g.tfname '_ALLCROSSFANGLE' ]); 
	eval(['load ' g.tffolder g.tfname '_times' ]); 
	eval(['load ' g.tffolder g.tfname '_freqs' ]); 
end;

% threshold activities (so that lines do not flash)
% -------------------------------------------------
try, 
	eval(['load ' g.tffolder g.tfname '_newERSP' ]);
	eval(['load ' g.tffolder g.tfname '_newITC' ]);
	eval(['load ' g.tffolder g.tfname '_newCROSSF' ]);
	eval(['load ' g.tffolder g.tfname '_newANGLE' ]);
catch,
	newERSP   = moviethresh( ALLERSP, 0.2, 4, 2);
	newITC    = moviethresh( ALLITC , 0.2, 4, 2);
	newCROSSF = moviethresh( ALLCROSSF, 0.2, 4, 2);
	newANGLE  = ALLCROSSFANGLE;
    %newANGLE  = revertangle2( ALLCROSSFANGLE, newCROSSF); % max angle
	
	eval(['save ' g.tffolder g.tfname '_newERSP   newERSP' ]);
	eval(['save ' g.tffolder g.tfname '_newITC    newITC' ]);
	eval(['save ' g.tffolder g.tfname '_newCROSSF newCROSSF' ]);
	eval(['save ' g.tffolder g.tfname '_newANGLE  newANGLE' ]);

	% this second call is innefective (but it is usefull for you to check that changes have been applied)
	%newERSP   = moviethresh( newERSP, 0.1, 4, 2);
	%newITC    = moviethresh( newITC , 0.1, 4, 2);
	%newCROSSF = moviethresh( newCROSSF, 0.1, 4, 2);
	%newANGLE  = revertangle2( newANGLE, newCROSSF);
end;

% frequencies to plot movie
% -------------------------
if isempty(g.freqs)
    g.freqindices = [1:length(freqs)];
else
    g.freqindices = [];
    for index = 1:length(g.freqs)
        [tmpfreq minfreq] = min(abs(freqs - g.freqs(index)));
        g.freqindices = [ g.freqindices minfreq];
        fprintf('Found closest frequency for %3.2f Hz: %3.2f Hz\n', g.freqs(index), freqs(g.freqindices(index)));
    end;
end;
        
% movie title defaults
% --------------------
for index = 1:nbconditions
    alltitles{index} = sprintf('Condition %d', index);
end;
if strcmpi(g.diffmovie, 'on')
    alltitles{3} = 'Cond1 - Cond2';
end

% pern10 movie parameters
% -----------------------
if isstr(g.movparams)& strcmpi(g.movparams, 'mriside')
        
    % -------------------
    % movie from the side
    % -------------------
    if isempty(g.coordinates)
        coordinates = founddipoles(ALLEEG, g.comps);
        [tmp plotorder] = sort( coordinates(:,1) );
        coordinates = coordinates(:, [2 3]); % remove X
    else
        plotorder   = g.showcomps;
        coordinates = g.coordinates;
    end;
    brainmovieoptions = { 'plotorder', plotorder(g.showcomps), ... 
                         'resolution', 'low', ...
                        'coordinates', coordinates, ...
                        'circfactor', g.circfactor, ...
                        'xlimaxes', [-1.15 1.15], ...
                        'ylimaxes', [-0.8 1.5], ...
                        'title', '', ...
                        'rthistloc', [9 9 1.3 1], ...
                        'envylabel', 'uV', ...
                        'visible', 'on', ...
                        'crossfphasespeed', 'off', ...
                        'head', 'mriside.pcx', ...
                        'crossfphaseunit', 'radian', ...
                        'size', [350 400], ...
                        'condtitleformat', { 'fontsize', 14, 'fontweight', 'bold'}, ...
                        'condtitle', alltitles };
elseif isstr(g.movparams) & strcmpi(g.movparams, 'mritop')
    
    % ------------------
    % movie from the top
    % ------------------
    if isempty(g.coordinates)
        coordinates = founddipoles(ALLEEG, g.comps);
        [tmp plotorder] = sort( coordinates(:,3) );
        coordinates = coordinates(:, [1 2]); % remove Z
    else
        plotorder   = g.showcomps;
        coordinates = g.coordinates;
    end;
    
    brainmovieoptions = {  'plotorder',  plotorder(g.showcomps), ...
                         'resolution', 'low', ...
                        'coordinates', coordinates, ...
                        'circfactor', g.circfactor, ...
                        'xlimaxes', [-1.1 1.1], ...
                        'ylimaxes', [-1.1 1.1], ...
                        'envylabel', 'uV', ...
                        'rthistloc', [9 9 1.3 1], ...
                        'visible', 'on', ...
                        'crossfphasespeed', 'off', ...
                        'head', '/data/common/matlab/toolbox2/mritop.pcx', ...
                        'crossfphaseunit', 'radian', ...
                        'size', [350 400], ...
                        'condtitleformat', { 'fontsize', 14, 'fontweight', 'bold'}, ...
                        'square', 'off', ...
                        'condtitle', alltitles };
elseif isstr(g.movparams) & strcmpi(g.movparams, 'mrirear')
    
    % ------------------
    % movie from the rear
    % ------------------
    if isempty(g.coordinates)
        coordinates = founddipoles(ALLEEG, g.comps);
        [tmp plotorder] = sort( coordinates(:,2) );
        coordinates = coordinates(:, [1 3]); % remove Z
    else
        plotorder   = g.showcomps;
        coordinates = g.coordinates;
    end;
    
    brainmovieoptions = {  'plotorder',  plotorder(g.showcomps), ...
                         'resolution', 'low', ...
                        'coordinates', coordinates, ...
                        'circfactor', g.circfactor, ...
                        'xlimaxes', [-1.1 1.1], ...
                        'ylimaxes', [-1.1 1.1], ...
                        'envylabel', 'uV', ...
                        'rthistloc', [9 9 1.3 1], ...
                        'visible', 'on', ...
                        'crossfphasespeed', 'off', ...
                        'head', '/data/common/matlab/toolbox2/mrirear.pcx', ...
                        'crossfphaseunit', 'radian', ...
                        'size', [350 400], ...
                        'condtitleformat', { 'fontsize', 14, 'fontweight', 'bold'}, ...
                        'square', 'off', ...
                        'condtitle', alltitles };
elseif isstr(g.movparams)
    error('Movparams template can only be ''mritop'' and ''mriside''');
else
    % ----------------------------------------------------------------
    % custom movie -> g.movparams contains cell array of movie options
    % ----------------------------------------------------------------
    brainmovieoptions = { 'condtitle' alltitles 'coordinates', g.coordinates, ...
                        'circfactor', g.circfactor, ...
                        g.movparams{:}};
end;
if strcmp(g.oneframe, 'on')
   brainmovieoptions = { brainmovieoptions{:} 'frames' [1] };
end;
if ~isempty(g.addmovparams)
    brainmovieoptions = { brainmovieoptions{:} g.addmovparams{:} };
end;

% data enveloppe
% --------------
for index = 1:nbconditions
    allenv(:,:,index) = env(mean(ALLEEG(index).data,3), [min(times) max(times)], times);
end;
if strcmpi(g.diffmovie, 'on')
	allenv(:,:,3) = env(mean(ALLEEG(1).data,3)-mean(ALLEEG(2).data,3), [min(times) max(times)], times);
end;
brainmovieoptions = { brainmovieoptions{:} 'envelope' allenv };

% plot polarity
% -------------
if strcmpi(g.diffmovie, 'on')
    brainmovieoptions = { brainmovieoptions{:} 'polarity' 'posneg' };
end;

% get reaction time
% -----------------
if ~isempty(g.eventprob)
    eventcellarray = {};
    for index = 1:nbconditions
        eventcellarray{index} = eeg_getepochevent(ALLEEG(index),g.eventprob );
    end;
	eventcellarray{nbconditions+1} = [];
    brainmovieoptions = { brainmovieoptions{:} 'rt' eventcellarray };
end;

% BRAINMOVIE 
% ----------
brainmovieoptions  =removedup(brainmovieoptions);
if ~strcmpi(g.mode, 'compute')
	origdir = pwd;
    
	for freqindex = g.freqindices
        [tmp1 tmp2] = mkdir('/', g.framefolder(2:end) );
        cd(g.framefolder);
        
		% Compute the MIN/MAX power 
		%--------------------------
        if ~strmatch('scalepower', brainmovieoptions(1:2:end)) % test if scalepower is defined
            tmpmin = 1000;
            tmpmax = -1000;
            for numcompo=1:length( g.comps ) 
                for condition=1:nbconditions 
                    tmpersp = ALLERSP{ numcompo, condition };
                    tmpmin  = min( tmpmin, min(	tmpersp(freqindex,:) ));
                    tmpmax  = max( tmpmax, max(	tmpersp(freqindex,:) ));
                end;
            end;
            brainmovieoptions = { brainmovieoptions{:} 'scalepower' [tmpmin tmpmax] };
        end;
        
        % Run brainmovie
        % --------------
        brainmovie( newERSP, newITC, newCROSSF, newANGLE, times, freqindex, g.showcomps, ...
                    brainmovieoptions{:}, 'framesout', fastif(strcmpi(g.quality, 'ultrafast'), 'ppm', 'fig'));  
        if strcmp(g.oneframe, 'on')
            disp('Only one frame generated');
            cd(origdir);
            return
        end;
        
        % Run makemovie
        % -------------
        if ~isempty(g.moviefolder)
             outname = sprintf('%s/%s%3.2f', g.moviefolder, g.moviename, freqs(freqindex));
        else outname = sprintf('%s%3.2f', g.moviename, freqs(freqindex));
        end;
        if strcmpi(g.quality, 'ultrafast')
            unix(sprintf('mkavi -file %s.avi image*.ppm', outname));
        else
            g.makemovie = removedup({ 'mode' g.quality g.makemovie{:} 'dir', g.framefolder, 'outname', outname });
            makemovie( { 'image' 1 length(times) 4 }, g.makemovie{:});
        end;
    end;    
	cd(origdir);
end;

return

% add a folder separator
% ----------------------
function str = addfinalsep(str)
    if isempty(str), return; end;
    if strcmpi(computer, 'PCWIN')
        if str(end) ~= '\', str(end+1) = '\'; end;
    elseif strcmpi(computer, 'MAC')
        if str(end) ~= ':', str(end+1) = ':'; end;
    else
        if str(end) ~= '/', str(end+1) = '/'; end;
    end;
    
% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique(cella(1:2:end));
    cella = cella(sort(union(indices*2-1, indices*2)));
    
% get dipoles location
% --------------------
function coordinates = founddipoles(ALLEEG, comps)
    if ~isfield(ALLEEG, 'sources')
        error('Field ''sources'' containing dipole location does not exist');
    end;
    
    % searching sources
    % -----------------
    indexeeg = find(~cellfun('isempty', { ALLEEG.sources }));
    indexeeg = indexeeg(1);
    fprintf('Found besa sources in dataset number %d\n', indexeeg);
    
    if ~isfield(ALLEEG(indexeeg).sources, 'X')
        fprintf('No 3-D coordinates found, running besaplot ...\n');
        ALLEEG(indexeeg).sources = besaplot(ALLEEG(indexeeg).sources);
        close;
    end;        
    
    % scanning components
    % -------------------
    for index = 1:length(comps)
        indexcomp = find(cell2mat({ALLEEG(indexeeg).sources.component}) == comps(index));
        if isempty(indexcomp)
            error(['Component ' int2str( comps(index) ) ' not found in the besa equivalent dipole strcuture']);
        end;
        if length(indexcomp) > 1
            error(['Warning: 2 equivalent dipoles found for component ' int2str( comps(index) ) ': only considering the first one']);
        end;            
        coordinates(index,1) = ALLEEG(indexeeg).sources(indexcomp(1)).X;
        coordinates(index,2) = ALLEEG(indexeeg).sources(indexcomp(1)).Y;
        coordinates(index,3) = ALLEEG(indexeeg).sources(indexcomp(1)).Z;
    end;
