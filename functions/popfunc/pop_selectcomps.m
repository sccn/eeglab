% pop_selectcomps() - Display components with button to visualize their
%                  properties and label them for rejection.
% Usage:
%       >> OUTEEG = pop_selectcomps( INEEG, compnum );
%
% Inputs:
%   INEEG    - Input dataset
%   compnum  - vector of component numbers
%
% Output:
%   OUTEEG - Output dataset with updated rejected components
%
% Note:
%   if the function POP_REJCOMP is ran prior to this function, some 
%   fields of the EEG datasets will be present and the current function 
%   will have some more button active to tune up the automatic rejection.   
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_prop(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_selectcomps( EEG, compnum, fig );

COLREJ = '[1 0.6 0.6]';
COLACC = '[0.75 1 0.75]';
PLOTPERFIG = 35;

com = '';
if nargin < 1
	help pop_selectcomps;
	return;
end;	

if nargin < 2
    uilist = { { 'style' 'text' 'string' 'Components to plot:' } ...
               { 'style' 'edit' 'string'  ['1:' int2str(size(EEG.icaweights,1)) ] } ...
               {} ...
               { 'style' 'text' 'string' [ 'Note: in the next interface, click on buttons to see' char(10) ... 
                                           'component properties and label them for rejection.' char(10) ...
                                           'To actually reject labelled components use menu item' char(10) ...
                                           '"Tools > Remove components" or use STUDY menus.' ] } };
                                       
    result = inputgui('uilist', uilist, 'geometry', { [1 1] 1 1 }, 'geomvert', [1 0.3 3], 'title', 'Reject comp. by map -- pop_selectcomps');
    if isempty(result), return; end
    compnum = eval( [ '[' result{1} ']' ]);

    if length(compnum) > PLOTPERFIG
        ButtonName=questdlg2(strvcat(['More than ' int2str(PLOTPERFIG) ' components so'],'this function will pop-up several windows'), ...
                             'Confirmation', 'Cancel', 'OK','OK');
        if ~isempty( strmatch(lower(ButtonName), 'cancel')), return; end
    end

end
fprintf('Drawing figure...\n');
currentfigtag = ['selcomp' num2str(rand)]; % generate a random figure tag

if length(compnum) > PLOTPERFIG
    for index = 1:PLOTPERFIG:length(compnum)
        pop_selectcomps(EEG, compnum([index:min(length(compnum),index+PLOTPERFIG-1)]));
    end

    com = [ 'pop_selectcomps(EEG, ' vararg2str(compnum) ');' ];
    return;
end

if isempty(EEG.reject.gcompreject)
	EEG.reject.gcompreject = zeros( size(EEG.icawinv,2));
end
try, icadefs; 
catch, 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end

% set up the figure
% -----------------
column =ceil(sqrt( length(compnum) ))+1;
rows = ceil(length(compnum)/column);
if ~exist('fig','var')
	figure('name', [ 'Reject components by map - pop_selectcomps() (dataset: ' EEG.setname ')'], 'tag', currentfigtag, ...
		   'numbertitle', 'off', 'color', BACKCOLOR);
	set(gcf,'MenuBar', 'none');
	pos = get(gcf,'Position');
	set(gcf,'Position', [pos(1) 20 800/7*column 600/5*rows*1.2]);
    incx = 120;
    incy = 110;
    sizewx = 100/column;
    if rows > 2
        sizewy = 90/rows;
	else 
        sizewy = 80/rows;
    end
    pos = get(gca,'position'); % plot relative to current axes
	hh = gca;
	q = [pos(1) pos(2) 0 0];
	s = [pos(3) pos(4) pos(3) pos(4)]./100;
	axis off;
end

% figure rows and columns
% -----------------------  
if EEG.nbchan > 64
    disp('More than 64 electrodes: electrode locations not shown');
    plotelec = 0;
else
    plotelec = 1;
end
count = 1;
for ri = compnum
	if exist('fig','var')
		button = findobj('parent', fig, 'tag', ['comp' num2str(ri)]);
		if isempty(button) 
			error( 'pop_selectcomps(): figure does not contain the component button');
		end;	
	else
		button = [];
	end;		
		 
	if isempty( button )
		% compute coordinates
		% -------------------
		X = mod(count-1, column)/column * incx-10;  
        	Y = (rows-floor((count-1)/column))/rows * incy - sizewy*1.3;  

		% plot the head
		% -------------
		if ~strcmp(get(gcf, 'tag'), currentfigtag);
		    figure(findobj('tag', currentfigtag));
		end
		ha = axes('Units','Normalized', 'Position',[X Y sizewx sizewy].*s+q);
		if plotelec
		    topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
			      'off', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
		else
		    topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
			      'off', 'electrodes','off', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
		end
		
		% labels
		% -------------
		if isfield(EEG.etc, 'ic_classification')
			classifiers = fieldnames(EEG.etc.ic_classification);
			if ~isempty(classifiers)
				if ~exist('classifier_name', 'var') || isempty(classifier_name)
					if any(strcmpi(classifiers, 'ICLabel'));
						classifier_name = 'ICLabel';
					else
						classifier_name = classifiers{1};
					end
				else
					classifier_name = classifiers{strcmpi(classifiers, classifier_name)};
				end
				if ri == compnum(1) && size(EEG.icawinv, 2) ...
						~= size(EEG.etc.ic_classification.(classifier_name).classifications, 1)
					warning(['The number of ICs do not match the number of IC classifications. This will result in incorrectly plotted labels. Please rerun ' classifier_name])
				end
				[prob, classind] = max(EEG.etc.ic_classification.(classifier_name).classifications(ri, :));
				t = title(sprintf('%s : %.1f%%', ...
					EEG.etc.ic_classification.(classifier_name).classes{classind}, ...
					prob*100));
				set(t, 'Position', get(t, 'Position') .* [1 -1.2 1])
			end
		end
		axis square;

		% plot the button
		% ---------------
         if ~strcmp(get(gcf, 'tag'), currentfigtag);
             figure(findobj('tag', currentfigtag));
         end
		button = uicontrol(gcf, 'Style', 'pushbutton', 'Units','Normalized', 'Position',...
                           [X Y+sizewy sizewx sizewy*0.18].*s+q, 'tag', ['comp' num2str(ri)]);
        command = sprintf('pop_prop( EEG, 0, %d, gcbo, { ''freqrange'', [1 50] });', ri);
		set( button, 'callback', command );
	end
	set( button, 'backgroundcolor', eval(fastif(EEG.reject.gcompreject(ri), COLREJ,COLACC)), 'string', int2str(ri)); 	
	drawnow;
	count = count +1;
end

% draw the bottom button
% ----------------------
if ~exist('fig','var')
    if ~strcmp(get(gcf, 'tag'), currentfigtag);
        figure(findobj('tag', currentfigtag));
    end
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Cancel', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[-10 -10  15 sizewy*0.25].*s+q, 'callback', 'close(gcf); fprintf(''Operation cancelled\n'')' );
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Set threhsolds', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[10 -10  15 sizewy*0.25].*s+q, 'callback', 'pop_icathresh(EEG); pop_selectcomps( EEG, gcbf);' );
	if isempty( EEG.stats.compenta	), set(hh, 'enable', 'off'); end;	
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'See comp. stats', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[30 -10  15 sizewy*0.25].*s+q, 'callback',  ' ' );
	if isempty( EEG.stats.compenta	), set(hh, 'enable', 'off'); end;	
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Help', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[70 -10  15 sizewy*0.25].*s+q, 'callback', 'pophelp(''pop_selectcomps'');' );
	command = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);''); close(gcf)';
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[90 -10  15 sizewy*0.25].*s+q, 'callback',  command);
			% sprintf(['eeg_global; if %d pop_rejepoch(%d, %d, find(EEG.reject.sigreject > 0), EEG.reject.elecreject, 0, 1);' ...
		    %		' end; pop_compproj(%d,%d,1); close(gcf); eeg_retrieve(%d); eeg_updatemenu; '], rejtrials, set_in, set_out, fastif(rejtrials, set_out, set_in), set_out, set_in));
end

com = [ 'pop_selectcomps(EEG, ' vararg2str(compnum) ');' ];
return;		
