% pop_selectcomps() - Display components with button to vizualize their
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
    promptstr = { 'Components to plot:' };
    initstr   = { [ '1:' int2str(size(EEG.icaweights,1)) ] };
    
    result = inputdlg2(promptstr, 'Reject comp. by map -- pop_selectcomps',1, initstr);
    if isempty(result), return; end;
    compnum = eval( [ '[' result{1} ']' ]);

    if length(compnum) > PLOTPERFIG
        ButtonName=questdlg2(strvcat(['More than ' int2str(PLOTPERFIG) ' components so'],'this function will pop-up several windows'), ...
                             'Confirmation', 'Cancel', 'OK','OK');
        if ~isempty( strmatch(lower(ButtonName), 'cancel')), return; end;
    end;

end;
fprintf('Drawing figure...\n');
currentfigtag = ['selcomp' num2str(rand)]; % generate a random figure tag

if length(compnum) > PLOTPERFIG
    for index = 1:PLOTPERFIG:length(compnum)
        pop_selectcomps(EEG, compnum([index:min(length(compnum),index+PLOTPERFIG-1)]));
    end;

    com = [ 'pop_selectcomps(' inputname(1) ', ' vararg2str(compnum) ');' ];
    return;
end;

if isempty(EEG.reject.gcompreject)
	EEG.reject.gcompreject = zeros( size(EEG.icawinv,2));
end;
try, icadefs; 
catch, 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end;

% set up the figure
% -----------------
column =ceil(sqrt( length(compnum) ))+1;
rows = ceil(length(compnum)/column);
if ~exist('fig')
	figure('name', [ 'Reject components by map - pop_selectcomps() (dataset: ' EEG.setname ')'], 'tag', currentfigtag, ...
		   'numbertitle', 'off', 'color', BACKCOLOR);
	set(gcf,'MenuBar', 'none');
	pos = get(gcf,'Position');
	set(gcf,'Position', [pos(1) 20 800/7*column 600/5*rows]);
    incx = 120;
    incy = 110;
    sizewx = 100/column;
    if rows > 2
        sizewy = 90/rows;
	else 
        sizewy = 80/rows;
    end;
    pos = get(gca,'position'); % plot relative to current axes
	hh = gca;
	q = [pos(1) pos(2) 0 0];
	s = [pos(3) pos(4) pos(3) pos(4)]./100;
	axis off;
end;

% figure rows and columns
% -----------------------  
if EEG.nbchan > 64
    disp('More than 64 electrodes: electrode locations not shown');
    plotelec = 0;
else
    plotelec = 1;
end;
count = 1;
for ri = compnum
	if exist('fig')
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
			disp('Aborting plot');
			return;
		end;
		ha = axes('Units','Normalized', 'Position',[X Y sizewx sizewy].*s+q);
        if plotelec
            topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
                      'off', 'style' , 'fill', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
        else
            topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
                      'off', 'style' , 'fill','electrodes','off', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
        end;
		axis square;

		% plot the button
		% ---------------
		button = uicontrol(gcf, 'Style', 'pushbutton', 'Units','Normalized', 'Position',...
                           [X Y+sizewy sizewx sizewy*0.25].*s+q, 'tag', ['comp' num2str(ri)]);
		command = sprintf('pop_prop( %s, 0, %d, %3.15f, { ''freqrange'', [1 50] });', inputname(1), ri, button);
		set( button, 'callback', command );
	end;
	set( button, 'backgroundcolor', eval(fastif(EEG.reject.gcompreject(ri), COLREJ,COLACC)), 'string', int2str(ri)); 	
	drawnow;
	count = count +1;
end;

% draw the bottom button
% ----------------------
if ~exist('fig')
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Cancel', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[-10 -10  15 sizewy*0.25].*s+q, 'callback', 'close(gcf); fprintf(''Operation cancelled\n'')' );
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Set threhsolds', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[10 -10  15 sizewy*0.25].*s+q, 'callback', 'pop_icathresh(EEG); pop_selectcomps( EEG, gcbf);' );
	if isempty( EEG.stats.compenta	), set(hh, 'enable', 'off'); end;	
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'See comp. stats', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[30 -10  15 sizewy*0.25].*s+q, 'callback',  ' ' );
	if isempty( EEG.stats.compenta	), set(hh, 'enable', 'off'); end;	
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'See projection', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[50 -10  15 sizewy*0.25].*s+q, 'callback', ' ', 'enable', 'off'  );
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Help', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[70 -10  15 sizewy*0.25].*s+q, 'callback', 'pophelp(''pop_selectcomps'');' );
	command = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);''); close(gcf)';
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[90 -10  15 sizewy*0.25].*s+q, 'callback',  command);
			% sprintf(['eeg_global; if %d pop_rejepoch(%d, %d, find(EEG.reject.sigreject > 0), EEG.reject.elecreject, 0, 1);' ...
		    %		' end; pop_compproj(%d,%d,1); close(gcf); eeg_retrieve(%d); eeg_updatemenu; '], rejtrials, set_in, set_out, fastif(rejtrials, set_out, set_in), set_out, set_in));
end;

com = [ 'pop_selectcomps(' inputname(1) ', ' vararg2str(compnum) ');' ];
return;		
