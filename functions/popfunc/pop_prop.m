% pop_prop() - plot the properties of a channel or of an independent
%              component. 
%
% Usage:
%   >> pop_prop( EEG, typecomp); % pops up a query window 
%   >> pop_prop( EEG, typecomp, chan, winhandle);
%
% Inputs:
%   EEG        - dataset structure (see EEGGLOBAL)
%   typecomp   - [0|1] compute electrode property (1) or component 
%                property (0). Default is 1.
%   chan       - channel or component number
%   winhandle  - if this parameter is present, buttons for the
%                rejection of the components are drawn. If 
%                non-zeros, this parameter is used to backpropagate
%                the color of the rejection button.
% 
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_runica(), eeglab()

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
% Revision 1.3  2002/09/04 23:46:30  arno
% debug
%
% Revision 1.2  2002/09/04 23:18:35  arno
% addding typecomp for output command
%
% Revision 1.1  2002/09/04 23:13:59  arno
% Initial revision
%
% Revision 1.17  2002/09/04 23:12:51  arno
% channel property
%
% Revision 1.16  2002/08/20 21:45:18  arno
% title/button overlap
%
% Revision 1.15  2002/08/12 22:41:03  arno
% [A[Adebug help button
%
% Revision 1.14  2002/08/12 21:01:15  arno
% text
%
% Revision 1.13  2002/08/12 16:24:15  arno
% [6~inputdlg2
%
% Revision 1.12  2002/08/12 14:59:42  arno
% button color
%
% Revision 1.11  2002/08/12 01:30:43  arno
% color
%
% Revision 1.10  2002/08/11 20:49:49  arno
% color
%
% Revision 1.9  2002/07/29 00:42:00  arno
% debugging
%
% Revision 1.8  2002/07/18 16:03:44  arno
% using spectopo to compute spectrum
%
% Revision 1.7  2002/07/12 16:46:19  arno
% accept button update
%
% Revision 1.6  2002/04/23 21:14:27  scott
% [same] -sm
%
% Revision 1.5  2002/04/23 21:13:02  scott
% test for EEG.icaact -> EEG.icaweights  -sm
%
% Revision 1.4  2002/04/07 20:24:47  scott
% adjusting erpimage title position -sm
%
% Revision 1.3  2002/04/07 20:01:44  scott
% dealing with title() in erpimage with 'erp' -sm
%
% Revision 1.2  2002/04/07 19:48:25  scott
% added 'erp' to erpimage -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision

% hidden parameter winhandle

% 01-25-02 reformated help & license -ad 
% 02-17-02 removed event index option -ad
% 03-17-02 debugging -ad & sm
% 03-18-02 text settings -ad & sm
% 03-18-02 added title -ad & sm

function com = pop_prop(EEG, typecomp, numcompo, winhandle)

com = '';
if nargin < 1
	help plotcomp;
	return;   
end;
if nargin == 1
	typecomp ==1;
end;
if typecomp == 0 & isempty(EEG.icaweights)
   error('No ICA weights recorded for this set, first run ICA');
end;   
if nargin == 2
	promptstr    = { fastif(typecomp,'Channel number to plot:','Component number to plot:') };
	inistr       = { '1' };
	result       = inputdlg2( promptstr, 'Component properties -- pop_prop()', 1,  inistr, 'pop_prop');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
   
	numcompo   = eval( result{1} );
end;

if numcompo < 1 | numcompo > EEG.nbchan
   error('Component index out of range');
end;   

% assumed input is numcompo
% -------------------------
try, icadefs; 
catch, 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end;
basename = [fastif(typecomp,'Channel ', 'Component ') int2str(numcompo) ];

fh = figure('name', [basename ' properties'], 'color', BACKCOLOR, 'numbertitle', 'off');
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2)-500+pos(4) 500 500]);
pos = get(gca,'position'); % plot relative to current axes
hh = gca;
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]./100;
axis off;

% plotting topoplot
% -----------------
h = axes('Units','Normalized', 'Position',[-10 60 40 42].*s+q);
%topoplot( EEG.icawinv(:,numcompo), EEG.chanlocs); axis square; 
if typecomp == 1
	topoplot( numcompo, EEG.chanlocs, 'electrodes','off', 'style', 'blank', 'emarkersize1chan', 10); axis square;
else
	topoplot( EEG.icawinv(:,numcompo), EEG.chanlocs, 'shading', 'interp', 'numcontour', 3); axis square;
end;
title([ basename fastif(typecomp, ' location', ' topography')], 'fontsize', 14); 

% plotting erpimage
% -----------------
h = axes('Units','Normalized', 'Position',[45 62 48 38].*s+q);
if EEG.trials > 1
	eeg_options; 
	title([basename ' activity'], 'fontsize', 14); % put title at top of erpimage
	axis off
	h = axes('Units','Normalized', 'Position',[45 62 48 38].*s+q);
	EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
	if typecomp == 1
		erpimage( EEG.data(numcompo,:), ones(1,EEG.trials)*10000, EEG.times , '', 3, 1, 'caxis', 1/4, 'cbar','erp');   
	else
		if option_computeica  
			erpimage( EEG.icaact(numcompo,:), ones(1,EEG.trials)*10000, EEG.times , '', 3, 1, 'caxis', 1/4, 'cbar','erp', 'yerplabel', '');   
		else
			icaacttmp = (EEG.icaweights(numcompo,:)*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
			erpimage( icaacttmp, ones(1,EEG.trials)*10000, EEG.times, '', 3, 1, 'caxis', 1/4, 'cbar','erp', 'yerplabel', '');   
		end;
	end;
else
	axis off;
	text(0.1, 0.3, [ 'No erpimage plotted' 10 'for continuous data']);
end;	

% plotting spectrum
% -----------------
if exist('winhandle')
	h = axes('units','normalized', 'position',[5 10 95 35].*s+q);
else
	h = axes('units','normalized', 'position',[5 0 95 40].*s+q);
end;
%h = axes('units','normalized', 'position',[45 5 60 40].*s+q);
try
	eeg_options; 
	if typecomp == 1
		spectopo( EEG.data(numcompo,:), EEG.pnts, EEG.srate );
	else 
		if option_computeica  
			spectopo( EEG.icaact(numcompo,:), EEG.pnts, EEG.srate );
		else
			if exist('icaacttmp')~=1, 
				icaacttmp = (EEG.icaweights(numcompo,:)*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts); 
			end;
			spectopo( icaacttmp, EEG.pnts, EEG.srate );
		end;
	end;
	set(gca, 'xlim', [0 min(50, EEG.srate/2)]);
	set( get(gca, 'ylabel'), 'string', 'Magnitude (dB)', 'fontsize', 14); 
	set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 14); 
	title('Activity power spectrum', 'fontsize', 14); 
catch
	axis off;
	text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
end;	
	
% display buttons
% ---------------

if exist('winhandle')
	COLREJ = '[1 0.6 0.6]';
	COLACC = '[0.75 1 0.75]';
	% CANCEL button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-10 -10 15 6].*s+q, 'callback', 'close(gcf);');

	% VALUE button
	% -------------
	hval  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [15 -10 15 6].*s+q);

	% REJECT button
	% -------------
	status = EEG.reject.gcompreject(numcompo);
	hr = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', eval(fastif(status,COLREJ,COLACC)), ...
				'string', fastif(status, 'REJECT', 'ACCEPT'), 'Units','Normalized', 'Position', [40 -10 15 6].*s+q, 'userdata', status, 'tag', 'rejstatus');
	command = [ 'set(gcbo, ''userdata'', ~get(gcbo, ''userdata''));' ...
				'if get(gcbo, ''userdata''),' ...
				'     set( gcbo, ''backgroundcolor'',' COLREJ ', ''string'', ''REJECT'');' ...
				'else ' ...
				'     set( gcbo, ''backgroundcolor'',' COLACC ', ''string'', ''ACCEPT'');' ...
				'end;' ];					
	set( hr, 'callback', command); 

	% HELP button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [65 -10 15 6].*s+q, 'callback', 'pophelp(''pop_prop'');');

	% OK button
	% ---------
 	command = [ 'global EEG;' ...
 				'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata'');' ...
 				'EEG.reject.gcompreject(' num2str(numcompo) ') = tmpstatus;' ]; 
	if winhandle ~= 0
	 	command = [ command ...
	 				sprintf('if tmpstatus set(%3.15f, ''backgroundcolor'', %s); else set(%3.15f, ''backgroundcolor'', %s); end;', ...
					winhandle, COLREJ, winhandle, COLACC)];
	end;				
	command = [ command 'close(gcf); clear tmpstatus' ];
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);

	% draw the figure for statistical values
	% --------------------------------------
	index = num2str( numcompo );
	command = [ ...
		'figure(''MenuBar'', ''none'', ''name'', ''Statistics of the component'', ''numbertitle'', ''off'');' ...
		'' ...
		'pos = get(gcf,''Position'');' ...
		'set(gcf,''Position'', [pos(1) pos(2) 340 340]);' ...
		'pos = get(gca,''position'');' ...
		'q = [pos(1) pos(2) 0 0];' ...
		's = [pos(3) pos(4) pos(3) pos(4)]./100;' ...
		'axis off;' ...
		''  ...
		'txt1 = sprintf(''(\n' ...
						'Entropy of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ' AND                 \t\t\t----\n\n' ...
					    'Kurtosis of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ') OR                  \t\t\t----\n\n' ...
					    'Kurtosis distibution \t\t\t%2.2f\n' ...
					    '> Rejection threhold\t\t\t%2.2f\n\n' ...
					    '\n' ...
					    'Current thesholds sujest to %s the component\n\n' ...
					    '(after manually accepting/rejecting the component, you may recalibrate thresholds for future automatic rejection on other datasets)'',' ...
						'EEG.stats.compenta(' index '), EEG.reject.threshentropy, EEG.stats.compkurta(' index '), ' ...
						'EEG.reject.threshkurtact, EEG.stats.compkurtdist(' index '), EEG.reject.threshkurtdist, fastif(EEG.reject.gcompreject(' index '), ''REJECT'', ''ACCEPT''));' ...
		'' ...				
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-11 4 117 100].*s+q, ''Style'', ''frame'' );' ...
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-5 5 100 95].*s+q, ''String'', txt1, ''Style'',''text'', ''HorizontalAlignment'', ''left'' );' ...
		'h = uicontrol(gcf, ''Style'', ''pushbutton'', ''string'', ''Close'', ''Units'',''Normalized'', ''Position'', [35 -10 25 10].*s+q, ''callback'', ''close(gcf);'');' ...
		'clear txt1 q s h pos;' ];
	set( hval, 'callback', command); 
	if isempty( EEG.stats.compenta )
		set(hval, 'enable', 'off');
	end;
	
	com = sprintf('pop_prop( %s, %d, %d, 0);', inputname(1), typecomp, numcompo );
else
	com = sprintf('pop_prop( %s, %d, %d);', inputname(1), typecomp, numcompo);
end;

return;

