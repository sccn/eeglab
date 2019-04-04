% rejstatepoch() - reject bad eeg trials based a statistical measure. Can 
%                 be applied either to the raw eeg data or the ICA 
%                 component activity. This is an interactive function.
%
% Usage:
%   >> [ Irej, Irejdetails, n, threshold, thresholdg] = ...
%                         rejstatepoch( signal, rej, 'key1', value1...);
%
% Inputs:
%   signal     - 3 dimensional data array channel x points x trials 
%                (instead of channels, one might also use independent
%                components).
%   rej        - rejection array, one value per trial and per channel or
%                component (size channel x trials). 
%                By default, values are normalized by this 
%                function and trehshold is expressed in term of standard 
%                deviation of the mean.
%
% Optional inputs:
%   'plot'       - ['on'|'off'] interactive mode or just rejection. 
%                  In the interactive mode, it plots the normalized 
%                  entropy of original signal (blue) and the limits 
%                  (red) (default:'on')            
%   'threshold'  - percentage error threshold (default 1-0.25/nb_trials)
%                  for individual trials of individual channel/component.
%                  This treshold is expressed in term of standart 
%                  deviation from the mean (default 5).
%   'global'     - ['on'|'off'], also perform threshold on the global
%                  measure (by default, the mean over all channel or
%                  electrodes).
%   'rejglob'    - rejection array, one value per trials. Use this 
%                  argument when the global measure for all channels or
%                  does not correspond to their mean.
%   'thresholdg' - global threshold for the reunion of all channels
%                  or components. By default, it is equal to 'threshold'.
%   'normalize'  - ['on'|'off'], normalize values before applying the 
%                  treshold. Default is 'on'.          
%   'plotcom'    - sting command to plot single trials. Default 
%                  is none.
%   'title'      - title of the graph. Default is none.
%   'labels'     - labels for electrodes (array not cell).
%
% Outputs:
%   Irej        - indexes of trials to be rejected
%   Irejdetails - array for rejected components or channel (nb_rejected x
%                 nb_channel or nb_rejected x nb_components)
%   n           - number of trials rejected
%   threshold   - percentage error threshold 
%   thresholdg  - percentage error threshold for global rejection 
%
% See also: eeglab()

% Algorithm:
%   normalise the measure given as input and reject trials based on 
%   an uniform distribution of the data

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

% userdata
% gcf : plotsig - signal to plot in the pop_out window
%       pnts - number of points per epoch
%       Irej - rejection of trials
%       Irejelec - rejection of trials and electrodes
%       labels - labels for curves 
% plotwin : rej1 - rejection array (with all electrodes x trials)
%           rej2 - global rejection array (size trials)
%           thr1 - threshold (electrodes)
%           thr2 - threshold global

function [ Irej, Irejdetails, n, threshold, thresholdg] = rejstatepoch( signal, ...
		rej, varargin); % pnts, th_E, th_rejg, command, commandplot, typetitle, E, rejg);

if nargin < 1
	help rejstatepoch;
	return;
end

if ~ischar( signal )

	if nargin < 2
		help rejstatepoch;
		return;
	end

	if ~isempty( varargin ), g=struct(varargin{:}); 
	else g= []; end
	try, g.plot; 			catch, g.plot='on'; end
	try, g.threshold; 		catch, g.threshold=5; end
	try, g.thresholdg; 		catch, g.thresholdg=5; end
	try, g.global; 			catch, g.global='on'; end
	try, g.rejglob; 		catch, g.rejglob=[]; end
	try, g.normalize; 		catch, g.normalize='on'; end
	try, g.plotcom; 		catch, g.plotcom=''; end
	try, g.title;	 		catch, g.title=''; end
	try, g.labels;	 		catch, g.labels=''; end

	g.rej = rej;
	clear rej
	switch lower(g.plot)
		case {'on', 'off'} ;  
		otherwise disp('Error: Plot must be either ''on'' or ''off'''); return;
	end;	
	switch lower(g.global)
		case {'on', 'off'} ;  
		otherwise disp('Error: Global must be either ''on'' or ''off'''); return;
	end;	
	switch lower(g.normalize)
		case {'on', 'off'} ;  
		otherwise disp('Error: Normalize must be either ''on'' or ''off'''); return;
	end;	
	if ~ischar(g.plotcom)
		disp('Error: Plotcom must be a string to evaluate'); return;
	end;	
	if ~ischar(g.title)
		disp('Error: Title must be a string'); return;
	end;	
	try, g.threshold*2;
		catch, disp('Error: Threhsold must be a number'); return;
	end;	
	try, g.thresholdg*2;
		catch, disp('Error: Threhsoldg must be a number'); return;
	end;	
	if length(g.threshold(:)) > 1
		disp('Error: Threhsold must be a single number'); return;
	end;	
	if length(g.thresholdg(:)) > 1
		disp('Error: Threhsoldg must be a single number'); return;
	end;	
	if ~isempty(g.rejglob)
		if length(g.rejglob) ~= size(g.rej,2)
			disp('Error: Rejglob must be have the same length as rej columns'); return;
		end
	else
		switch lower(g.global), case 'on', g.rejglob = sum(g.rej,1); end
	end;		
	if size(signal,3) ~= size(g.rej,2)
		disp('Error: Signal must be have the same number of element in 3rd dimension as rej have columns'); return;
	end
	if isempty(g.labels)
		for index = 1:size(g.rej,1)
			g.labels(index,:) = sprintf('%3d', index);
		end
		if ~isempty(g.rejglob)
			g.labels(index+2,:) = 'g. ';
		end
	end;		
			
	switch lower(g.normalize),  
		case 'on', 
			g.rej = (g.rej-mean(g.rej,2)*ones(1, size(g.rej,2)))./ (std(g.rej, 0, 2)*ones(1, size(g.rej,2)));
			switch lower(g.global), case 'on', g.rejglob = (g.rejglob(:)-mean(g.rejglob(:)))./ std(g.rejglob(:)); end
	end;	
	switch lower(g.global), case 'off',g.rejglob = []; end

	% plot the buttons
	% ----------------
	try, icadefs; catch, GUIBUTTONCOLOR = [0.8 0.8 0.8]; BACKCOLOR = [0.8 0.8 0.8]; end; 
	figure('color', BACKCOLOR);
	set(gcf, 'name', 'Rejectrials');
	pos = get(gca,'position'); % plot relative to current axes
	set( gca, 'tag', 'mainaxis');
	q = [pos(1) pos(2) 0 0];
	s = [pos(3) pos(4) pos(3) pos(4)]./100; % allow to use normalized position [0 100] for x and y
	axis('off');
	
	plotsig = sum(abs(signal(:,:)),1);
	set(gcf, 'userdata', { plotsig, size(signal, 2), [], [] }); % the two last arguments are the rejection

	% Create axis
	% -----------
	h6 = axes('Units','Normalized', 'tag',  'Plotwin', 'Position',[-10 12 120 84].*s+q);
	title(g.title);
	set( h6, 'userdata', { g.rej g.rejglob g.threshold g.thresholdg g.labels }); % g.rej was put twice because it is used to compute the global entropy

	% CANCEL button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Cancel', 'Units','Normalized','Position',...
			[-15 -6 15 6].*s+q, 'callback', 'close(gcf);', 'backgroundcolor', GUIBUTTONCOLOR);

	% Entropy component text and threshold
	% ------------------------------------
    makebutton( 'Single-channel', 'Plotwin'  ,  [3 -6 27 6].*s+q, [30 -6 15 6].*s+q, 3, g.threshold, GUIBUTTONCOLOR); 
    makebutton( 'All-channel', 'Plotwin'     ,  [50 -6  27 6].*s+q, [77 -6  15 6].*s+q, 4, g.thresholdg, GUIBUTTONCOLOR); 

	% EEGPLOT button
	% --------------
	if ~isempty(g.plotcom)
			h  = uicontrol(gcf, 'backgroundcolor', GUIBUTTONCOLOR, 'Style', 'pushbutton', 'string', 'EEGPLOT', 'Units','Normalized','Position',[95 -2 15 6].*s+q, ...
					'callback',['TMPEEG = get(gcbf, ''userdata'');' ...
								'TMPREJ = TMPEEG{3};' ...
								'TMPREJELEC = TMPEEG{4};' ...
								g.plotcom ] );
			posaccept = [95 -10 15 6];
	else	posaccept = [95 -6 15 6];
	end;								

	% ACCEPT button
	% -------------
	command = 'fprintf(''Rejection indices has been put in the matrix TMPREJ\n'')'; 
	haccept  = uicontrol(gcf, 'backgroundcolor', GUIBUTTONCOLOR, 'Style', 'pushbutton', 'string', 'UPDATE', 'Units','Normalized','Position', posaccept.*s+q, ...
					'callback', [	'set(gcbo, ''userdata'', 1);' ... %signal to signify termination 
									'TMPEEG = get(gcbf, ''userdata'');' ...
									'TMPREJ = TMPEEG{3};' ...
									'TMPREJELEC = TMPEEG{4};' ...
									command ] );

	command = [ 'entwin = get(gcbf, ''currentobject'');' ... 
				'tmptype = get(entwin,''type'');' ...
				'if tmptype(1) == ''l'' entwin = get(entwin, ''parent''); end;' ... 
				'tagwin = get(entwin,''tag'');' ... % either entropy or kurtosis 
				'switch tagwin,' ...
				' case ''Plotwin'',' ... % check if the user clicked on the right window
				'   alldata = get(gcf, ''userdata'');' ...
				'   plotsig = alldata{1};' ... 
				'   pnts = alldata{2};' ... 
				'   fig = figure(''position'', [100 300 600 400],''color'', [1 1 1]);' ...
				'   I = alldata{3};' ... 
				'   sweeps = size(plotsig,2) /pnts;' ... 
				'   h1 = axes(''parent'', fig, ''Units'',''Normalized'', ''Position'',[0.6 0.11 0.27 0.815]);' ...
				'   pos = get(entwin, ''currentpoint'');' ...
				'   component = round(pos(1) / sweeps + 0.5);' ... % determine the component number
				'   alldata = get(entwin, ''userdata'');' ...
		    	'   rej = alldata{1};' ... 
				'   if component <= size(rej,1)' ... % component 
				'        rej_threshold = alldata{3};' ... 
				'        component = max(1,min(component, size(rej,1))); ' ... 
		    	'   	 rej = rej(component, :);' ... 
				'        titlegraph = sprintf(''' g.title ' #%d'', component);' ... 
				'        colorgraph = ''b'';' ...
				'   else' ...                        % global 
		    	'   	 rej = alldata{2};' ... 
				'        rej_threshold = alldata{4};' ... 
				'        titlegraph = sprintf(''' g.title ' global'');' ... 
				'        colorgraph = ''g'';' ...
				'   end;' ...
				'   plot([1:length(rej)], rej, ''color'', colorgraph);' ... 
				'   ss = get(h1, ''xlim'');' ...
				'   set(h1, ''view'', [90 90]);' ...
				'   set(h1, ''xdir'', ''reverse'');' ...
				'   set(h1, ''XLim'', ss);' ...
				'   hold on;' ...  % plot component
				'   yl = get(h1, ''ylim'');' ...
				'   set(h1, ''xtickmode'', ''manual'', ''xtick'', sweeps/2, ''xticklabel'', component, ''xlim'', [ 0 sweeps ]);' ...
				'   title( titlegraph );' ...
				'   plot( get(h1, ''xlim''), [rej_threshold rej_threshold], ''r'');' ... % plot limit		  
				'   plot( get(h1, ''xlim''), [-rej_threshold -rej_threshold], ''r'');' ... % plot limit		  
				'   set(h1, ''xticklabel'', []);' ...
				'   hold on;' ...
				'   h2 = axes(''parent'', fig,''Units'',''Normalized'', ''Position'',[0.13 0.11 0.27 0.815]);' ...
				'   erpimage( plotsig, ones(1,size(plotsig,2)/pnts), [0:1/(pnts-1):1], '''', 3, 1, ''nosort'', ''noplot'');' ...
				'   title(''Currentset all chans''); xlabel(''''); ylabel(''''); ' ...
				'   set(gca, ''xticklabel'', []);' ...
				'   hold on;' ...
				'   h3 = axes(''parent'', fig,''Units'',''Normalized'', ''Position'',[0.45 0.11 0.1 0.815]);' ...
				'   if any(I ~= 0)' ...
				'      rejImage = (I'' * ones(1, 10))'';' ...
				'      imagesc( rejImage'' );' ...
				'      set(gca, ''ydir'', ''normal'');' ...
				'   end;' ...
				'   title(''Rejected (all)''); xlabel(''''); ylabel('''');' ... 
				'   set(gca, ''xticklabel'', [], ''yticklabel'', []);' ... 
				'end;' ...
				'clear fig tmptype tagwin alldata rej rejImage plotsig sweeps pnts rej_threshold ss q s h1 h2 h3 pos component yl;' ];

%				'      erpimage( rejImage(:)'', ones(1,size(I,2)), [0:1/(10-1):1], '''', 1, 0, ''nosort'', ''noplot'');' ...
	set(gcf, 'WindowButtonDownFcn', command);			

	rejstatepoch('draw');
	switch g.plot,
		case 'on', waitfor( haccept, 'userdata'); drawnow;
	end

	threshold  = g.threshold;
	thresholdg = g.thresholdg;
	Irej = [];
	Irejdetails = [];
	n = 0;
	try
		TMPEEG = get(gcf, 'userdata');
		Irej = TMPEEG{3};
		Irejdetails = TMPEEG{4};
		n = length(find(Irej == 1));

		plothandler = findobj( 'parent', gcf, 'tag', 'Plotwin');
		TMPEEG = get(plothandler, 'userdata');
		threshold = TMPEEG{3};
		thresholdg = TMPEEG{4};
		close(gcf);
	catch, end
else %if signal is a string draw everything

	% retreive data
	% -------------
	gcfdata = get(gcf, 'userdata');
	plotsig = gcfdata {1};
	pnts    = gcfdata {2};
    sweeps  = size(plotsig,2)/pnts;

	h6 = findobj('parent', gcf, 'tag', 'Plotwin');
	alldata = get(h6, 'userdata');
	g.rej       = alldata {1};
	g.rejg      = alldata {2};
	g.threshold   = alldata {3};
	g.thresholdg  = alldata {4};
	set(h6, 'userdata', alldata);

	nbchans = size(g.rej,1);

	% reject trials
	% -------------
	rejelec = abs(g.rej) > g.threshold;
	rej  = max(rejelec,[],1);
	n1 = sum(rej(:));
	if ~isempty( g.rejg )
		rej2 = abs(g.rejg) > g.thresholdg;
		n2 = sum(rej2(:));
		rej = rej | rej2(:)';
	end
	fprintf('%d trials rejected (single:%d, all:%d)\n', sum(rej(:)), n1, n2);
	gcfdata {3} = rej;
	gcfdata {4} = rejelec;
	set(gcf, 'userdata', gcfdata);
	
	% plot the sorted entropy curve
	% -----------------------------
	plotstat( 'Plotwin');

end
return;

function plotstat( id );

	% plot the sorted entropy curve
	% -----------------------------
	h6 = findobj('parent', gcf, 'tag', id);
	axes(h6); cla;
	ttmp = get(gca, 'title');
	oldtitle = get(ttmp, 'string');

	% get datas
    % ---------
   	alldata = get(gca, 'userdata');
	g.rej       = alldata {1};
	g.rejg      = alldata {2};
	g.threshold   = alldata {3};
	g.thresholdg  = alldata {4};
	g.labels      = alldata {5};
	nbchans = size(g.rej,1);
	sweeps  = size(g.rej,2);

	% plot datas
    % ----------
	g.rej = g.rej'; plot(g.rej(:)); g.rej = g.rej'; 
	hold on;
	yl = get(gca, 'ylim');

	% plot vertival bars to separate components and the trehsold
	% ----------------------------------------------------------
	set( gca, 'tag',  id, 'ylimmode', 'manual');
	set(gca, 'xtickmode', 'manual', 'xtick', [0:sweeps:(size(g.rej(:),1)-1+2*sweeps)] + sweeps/2, ...
			 'xticklabel', g.labels, 'xlim', [ 0 (size(g.rej(:),1)-1+2*sweeps)]);
	plot( [1 size(g.rej(:),1)], [-g.threshold -g.threshold], 'r');	% plot threshold	  
	plot( [1 size(g.rej(:),1)], [g.threshold g.threshold], 'r');	% plot threshold	  

	if ~isempty(g.rejg) % plot global ?	 
		plot([size(g.rej(:),1)+sweeps:size(g.rej(:),1)+2*sweeps-1],  g.rejg(:), 'g');
		pp = patch([size(g.rej(:),1) size(g.rej(:),1) size(g.rej(:),1)+sweeps size(g.rej(:),1)+sweeps], [yl(1)-1 yl(2)+1 yl(2)+1 yl(1)-1], get(gcf, 'color'), 'clipping', 'off');
		set(pp, 'EdgeColor',  get(gcf, 'color'));
		plot( [size(g.rej(:),1)+sweeps length(g.rejg)+size(g.rej(:),1)+sweeps], [-g.thresholdg -g.thresholdg], 'r');	% plot threshold	  
		plot( [size(g.rej(:),1)+sweeps length(g.rejg)+size(g.rej(:),1)+sweeps], [g.thresholdg g.thresholdg], 'r');	% plot threshold	  
		plot([size(g.rej(:),1)+sweeps size(g.rej(:),1)+sweeps], yl, 'k');
	else
		pp = patch([size(g.rej(:),1) size(g.rej(:),1) size(g.rej(:),1)+2*sweeps size(g.rej(:),1)+2*sweeps], [yl(1)-1 yl(2)+1 yl(2)+1 yl(1)-1], get(gcf, 'color'), 'clipping', 'off');
		set(pp, 'EdgeColor',  get(gcf, 'color'));
	end
	for index = 0:sweeps:size(g.rej(:),1); 
		plot([index index], yl, 'k');
	end

	% restore properties
	title(oldtitle);
	set(h6, 'userdata', alldata);

return;

function makebutton( string, ax, pos1, pos2, userindex, init, GUIBUTTONCOLOR );
	h  = uicontrol(gcf , 'backgroundcolor', GUIBUTTONCOLOR, 'Style', 'radiobutton', 'string', string, 'value', fastif(init == 0, 0, 1), 'Units','Normalized', 'Position', pos1, ...
					'callback', [ 'textresh = findobj(''parent'', gcbf, ''tag'', ''' string ''');' ...
								  'checkstatus = get(gcbo, ''value'');' ...
								  'ax = findobj(''parent'', gcbf, ''tag'', ''' ax ''');' ...
								  'tmpdat = get(ax, ''userdata'');' ...
								  'if checkstatus' ... % change status of the textbox 
								  '   set(textresh, ''enable'', ''on'');' ...
								  '   tmpdat{' int2str(userindex) '} = str2num(get(textresh, ''string''));' ...
								  'else' ... 		
								  '   set(textresh, ''enable'', ''off'');' ...
								  '   tmpdat{' int2str(userindex) '} = 0;' ...
								  'end;' ...
								  'set(ax, ''userdata'' , tmpdat);' ...
								  'rejstatepoch(''draw'');' ...
								  'clear tmpdat ax checkstatus textresh;'  ] );
	h  = uicontrol(gcf, 'Style', 'edit', 'backgroundcolor', [1 1 1], 'tag', string, 'string', num2str(init), 'enable', fastif(init == 0, 'off', 'on'), 'Units','Normalized', 'Position', pos2, ...
					'callback', [ 	'ax = findobj(''parent'', gcbf, ''tag'', ''' ax ''');' ...
								    'tmpdat = get(ax, ''userdata'');' ...
									'tmpdat{' int2str(userindex) '} = str2num(get(gcbo, ''string''));' ...
									'set(ax, ''userdata'' , tmpdat);' ...
								    'rejstatepoch(''draw'');' ...
									'clear tmpdat ax;' ] );
return;
