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
%                structural side MRI image as background. Default is 'mriside'.
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
                            'showcomps'   'integer'       []                                       [];
                            'makemovie'   'cell'          {}                                       {};
                            'eventprob'   ''              []                                       [] });
if isstr(g), error(g); end;
                    
% checking parameters
% -------------------
if strcmpi(g.diffmovie, 'on') &	length(ALLEEG) ~= 2
    error('For difference movies: nedd exactly to process 2 datasets');
end;
if isempty(g.showcomps), g.showcomps = g.comps; end;
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
	newANGLE  = revertangle2( ALLCROSSFANGLE, newCROSSF); % max angle
	
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
        fprintf('Found closest frequency for %3.2f Hz: %3.2f Hz\n', g.freqs(index), tmpfreq);
    end;
end;
        
% pern10 movie parameters
% -----------------------
if strcmpi(g.mode, 'mriside')
    
    % -------------------
    % movie from the side
    % -------------------
    for index = 1:nbconditions
        alltittles{index} = sprintf('Condition %d', index);
    end;
    alltittles = strvcat(alltitles{:});
    if isempty(g.coordinates)
        if ~isfield(ALLEEG, 'sources')
            error('Field ''sources'' containing dipole location does not exist');
        end;
        for index = g.comps
            coordinates(index,1) = EEG.sources(index).Y;
            coordinates(index,2) = EEG.sources(index).Z;
        end;
    else
        coordinates = g.coordinates;
    end;
    brainmovieoptions = { 'resolution', 'low', ...
                        'coordinates', coordinates, ...
                        'circfactor', g.circfactor, ...
                        'xlimaxes', [-1.1 1.1], ...
                        'ylimaxes', [-1.1 1.1], ...
                        'title', '', ...
                        'rthistloc', [9 9 1.3 1], ...
                        'envelope', allenv, ...
                        'envylabel', 'uV', ...
                        'visible', 'on', ...
                        'crossfphasespeed', 'off', ...
                        'head', 'mrirot.pcx', ...
                        'crossfphaseunit', 'radian', ...
                        'size', [350 400], ...
                        'condtitleformat', { 'fontsize', 14, 'fontweight', 'bold'}, ...
                        'condtitle', alltittles(1:nbconditions,:) };
elseif strcmpi(g.mode, 'mritop')
    
    % ------------------
    % movie from the top
    % ------------------
    for index = 1:nbconditions
        alltittles{index} = sprintf('Condition %d', index);
    end;
    alltittles = strvcat(alltitles{:});
    
    if isempty(g.coordinates)
        if ~isfield(ALLEEG, 'sources')
            error('Field ''sources'' containing dipole location does not exist');
        end;
        for index = g.comps
            coordinates(index,1) = EEG.sources(index).X;
            coordinates(index,2) = EEG.sources(index).Y;
        end;
    else
        coordinates = g.coordinates;
    end;
    
    brainmovieoptions = { 'resolution', 'low', ...
                        'coordinates', coordinates, ...
                        'circfactor', g.circfactor, ...
                        'xlimaxes', [-1.1 1.1], ...
                        'ylimaxes', [-1.1 1.1], ...
                        'envelope', allenv, ...
                        'envylabel', 'uV', ...
                        'rthistloc', [9 9 1.3 1], ...
                        'visible', 'on', ...
                        'crossfphasespeed', 'off', ...
                        'head', '/data/common/matlab/toolbox2/mritop.pcx', ...
                        'crossfphaseunit', 'radian', ...
                        'size', [350 400], ...
                        'condtitleformat', { 'fontsize', 14, 'fontweight', 'bold'}, ...
                        'square', 'off', ...
                        'condtitle', alltittles(1:nbconditions,:) };
else 
    % ------------
    % custom movie
    % ------------
    brainmovieoptions = { 'coordinates', g.coordinates, ...
                        'circfactor', g.circfactor, ...
                        g.movparams{:}};
end;
if ~isempty(g.addmovparams)
    brainmovieoptions = { brainmovieoptions{:} g.addmovparams };
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
%if strcmp(g.oneframe, 'on')
%   brainmovieoptions = { 'frames' [1 2] };
%end;

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
		tmpmin = 1000;
		tmpmax = -1000;
		for numcompo=1:length( g.comps ) 
			for condition=1:nbconditions 
				tmpersp = ALLERSP{ numcompo, condition };
				tmpmin  = min( tmpmin, min(	tmpersp(freqindex,:) ));
				tmpmax  = max( tmpmax, max(	tmpersp(freqindex,:) ));
			end;
		end;
		
        % Run brainmovie
        % --------------
		brainmovie( newERSP, newITC, newCROSSF, newANGLE, times, freqindex, g.showcomps, ...
                    brainmovieoptions{:}, 'framesout', 'fig', 'scalepower', [tmpmin tmpmax] );  
        %if strcmp(g.oneframe, 'on')
        %    disp('Only one frame generated');
        %    return
        %end;
        
        % Run makemovie
        % -------------
        g.framefolder = pwd;
        if ~isempty(g.moviefolder)
            cd(g.moviefolder)
        else 
            cd(origdir);
        end;
        if length(g.freqindices) > 1, outname = g.moviename;
        else                          outname = sprintf('%s3.2f', g.moviename, freqs(g.freqindices(freqindex)));
        end;
        g.makemovie = removedup({ 'mode' 'fast' g.makemovie{:} 'dir', g.framefolder, 'outname', outname });
        makemovie( { 'image' 1 length(times) 4 }, g.makemovie{:});
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