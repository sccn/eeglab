% pop_loadbci() - import a BCI2000 ascii file into EEGLAB
%
% Usage:
%   >> OUTEEG = pop_loadbci( filename, srate );
%
% Inputs:
%   filename       - file name
%   srate          - sampling rate
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 9 July 2002
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.9  2003/04/10 18:05:22  arno
% default argument
%
% Revision 1.8  2002/11/13 17:09:12  scott
% help msg
%
% Revision 1.7  2002/10/15 17:00:49  arno
% drawnow
%
% Revision 1.6  2002/08/12 16:33:37  arno
% inputdlg2
%
% Revision 1.5  2002/07/11 17:18:15  arno
% updating header
%
% Revision 1.4  2002/07/10 23:44:17  arno
% debugging events .
%
% Revision 1.3  2002/07/09 23:30:50  arno
% selecting interesting events
%
% Revision 1.2  2002/07/08 19:28:15  arno
% changing function filter
%
% Revision 1.1  2002/07/08 19:20:53  arno
% Initial revision
%

function [EEG, command] = pop_loadbci(filename, srate); 
EEG = [];
command = '';

if nargin < 1
	% ask user
	[filename, filepath] = uigetfile('*.*', 'Choose a BCI file -- pop_loadbci'); 
    drawnow;
	if filename == 0 return; end;
	filename = [filepath '/' filename];
	promptstr    = { 'Sampling rate' };
	inistr       = { '256' };
	result       = inputdlg2( promptstr, 'Import BCI2000 data -- pop_loadbci()', 1,  inistr, 'pop_loadbci');
	if length(result) == 0 return; end;
	srate   = eval( result{1} );
end;

% try to read as matlab
EEG = eeg_emptyset;
fprintf('Pop_loadbci: importing BCI file...\n');
try
	load( filename, '-mat');
catch
	
	fid = fopen(filename, 'r');
	allcollumns = fgetl(fid);
	colnames = {};
	while ~isempty(deblank(allcollumns))
		[colnames{end+1} allcollumns] = strtok(allcollumns);
	end;
	tmpdata = fscanf(fid, '%f', Inf);
	tmpdata = reshape(tmpdata, length(colnames), length(tmpdata)/length(colnames));
	%tmpdata = tmpdata';
	
	% find electrode indices
	% ----------------------
	indices = [];
	for index = 1:length(colnames)
		if strcmp( colnames{index}(1:2), 'ch')
			indices = [indices index ];
		end;
	end;
	eventindices = setdiff(1:length(colnames), indices);
    
    % ask for which event to import
    % -----------------------------
    geom = {[0.7 0.7 0.7 0.7]};
    uilist = { { 'style' 'text' 'string' 'State name' 'fontweight' 'bold' } ...
               { 'style' 'text' 'string' 'Import state' 'fontweight' 'bold'  } ...
               { 'style' 'text' 'string' 'Adjust time'  'fontweight' 'bold' } ...
               { 'style' 'text' 'string' 'Attribute of'  'fontweight' 'bold' } };
    allfields = strvcat('none', colnames{eventindices});
    for index = 1:length(eventindices)
        if strcmpi( colnames{eventindices(index)}, 'sourcetime')
            geom   = { geom{:} [1.3 0.3 0.3 0.5 0.1 1] };
            uilist = { uilist{:} { 'style' 'text' 'string' colnames{eventindices(index)} } ...
                                 { } { } { 'style' 'text' 'string' 'ref' } {} ...
                                 { } };
        else 
            geom   = { geom{:} [1.3 0.3 0.3 0.3 0.3 1] };
            uilist = { uilist{:} { 'style' 'text' 'string' colnames{eventindices(index)} } };
            if ~isempty(findstr( lower(colnames{eventindices(index)}), 'time')) 
                uilist = { uilist{:} { 'style' 'checkbox' 'value' 1 } { } { 'style' 'checkbox' 'value' 1 } {} ...
                           { 'style' 'listbox' 'string' '' } };
            elseif ~isempty(findstr( lower(colnames{eventindices(index)}), 'type'))
                uilist = { uilist{:} { 'style' 'checkbox' 'value' 1 } { } { 'style' 'checkbox' } {}  ...
                         { 'style' 'listbox' 'string' allfields 'value' min(length(eventindices), 6) ...
                         'listboxtop' min(length(eventindices), 6) }};
            else
                uilist = { uilist{:} { 'style' 'checkbox' } { } { 'style' 'checkbox' } {} ...
                         { 'style' 'listbox' 'string' allfields } };
            end;
        end;
    end;
    result = inputgui( geom, uilist, 'pophelp(''pop_loadbci'')', 'Import BCI2000 data files - pop_loadbci()');
    if isempty(result), return; end;
    
    % decoding result
    % ---------------
    listimport = {};
    s = 0;
    for index = 0:length(eventindices)-1
        tmplist = {};
        if strcmpi( colnames{eventindices(index+1)}, 'SourceTime')
            s = 1;
        else
            i = index-s;
            if result{3*i+1}, tmplist{1} = colnames{eventindices(index+1)}; else tmplist{1} = []; end;
            if result{3*i+2}, tmplist{2} = 1; else tmplist{2} = 0; end;
            if result{3*i+3}~=1 && result{3*i+3}~=0, tmplist{3} = colnames{ eventindices(result{3*i+3}-1+s) }; end;
            if ~isempty(tmplist{1}), listimport = { listimport{:} 'event' tmplist}; end;
        end;
    end;
    
    % find indices
    % ------------
    indexsource =  strmatch('SourceTime', colnames);
    count = 1;
    for index = 2:2:length(listimport)
        tmpindmatch = strmatch(listimport{index}{1}, colnames);
        if ~isempty(tmpindmatch), indeximport(count) = tmpindmatch; 
        else error(['State ''' listimport{index}{1} ''' not found']); 
        end;
        if length( listimport{index} ) > 1
             adjust(count) = listimport{index}{2};
        else adjust(count) = 0;
        end;
        if length( listimport{index} ) > 2
            tmpindmatch = strmatch(listimport{index}{3}, colnames);
            if ~isempty(tmpindmatch), corresp(count) = tmpindmatch; 
            else error(['State ''' listimport{index}{3} ''' not found']); 
            end;
        else
            corresp(count) = 0;
        end;
        count = count+1;
    end;
    indeximport
    corresp
    adjust
    
    % find block size
    % ---------------
    tmpevent = find( diff(tmpdata(indexsource, :)) ~= 0);
    diffevent = tmpevent(2:end)-tmpevent(1:end-1);
    blocksize = unique(diffevent);
    if length(blocksize) > 1, error('Error in determining block size'); 
    else                      fprintf('Blocksize: %d\n', blocksize); end;
    
    % segregate type and latency
    % --------------------------
    tmpcorresp = find(corresp);
    if length(tmpcorresp) ~= length(intersect(corresp, indeximport))
        disp('Warning: correspondance problem, some information will be lost');
    end;
    indexcorresp =  indeximport(tmpcorresp);
    indeximport(tmpcorresp) = [];
    adjust(tmpcorresp)      = [];
    
	% process events
	% --------------
	fprintf('Pop_loadbci: importing events...\n');
	counte = 1; % event counter
	events(10000).latency = 0;
	for index = 1:length(indeximport)
		tmpevent = find( diff(tmpdata(indeximport(index), :)) ~= 0);
		tmpevent = tmpevent+1;
        tmpcorresp = find(indexcorresp == indeximport(index));
        if adjust(index)
             fprintf('Latency of event ''%s'' adjusted\n', colnames{ indeximport(index) });
        else fprintf('WARNING: Latency of event ''%s'' not adjusted (latency uncertainty %2.1f ms)\n', ...
                     colnames{ indeximport(index) }, blocksize/srate*1000);
        end;
        
		for tmpi = tmpevent
            curlatency  = tmpdata(indeximport(index), tmpi);
			if curlatency % non zero
                if ~isempty(tmpcorresp)
                    events(counte).type = colnames{ indexcorresp(tmpcorresp) };
                else 
                    events(counte).type = colnames{ indeximport(index) };
                end;
                if adjust(index)
                    baselatency = tmpdata(indexsource, tmpi); % note that this is the first bin a block
                    realtmpi    = tmpi+blocksize;             % jump to the end of the block+1
                    if curlatency < baselatency, curlatency = curlatency+65536; end; % in ms
                    events(counte).latency = realtmpi+(curlatency-baselatency)/1000*srate;
                    % there is still a potentially large error between baselatency <-> realtmpi
                else
                    events(counte).latency = tmpi+(blocksize-1)/2;
                end;
				counte = counte+1;
			end;
		end;
	end;
	
	EEG.data = tmpdata(indices,:);
	EEG.nbchan = size(EEG.data, 1);
	EEG.srate  = srate;
	EEG = eeg_checkset(EEG);
	EEG.event = events(1:counte-1);	
	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
	
	% add up or down events
	% ---------------------
% $$$     for index=1:length(EEG.event)
% $$$ 		if strcmp(EEG.event(index).type(1:6), 'Target')
% $$$ 			targetcode = str2num(EEG.event(index).type(end));
% $$$ 			if targetcode == 1
% $$$ 				EEG.event(index).type = 'toptarget';
% $$$ 			else
% $$$ 				EEG.event(index).type = 'bottomtarget';
% $$$ 			end;
% $$$ 		else 
% $$$ 			if strcmp(EEG.event(index).type(1:6), 'Result')
% $$$ 				resultcode = str2num(EEG.event(index).type(end));
% $$$ 				if resultcode == 1
% $$$ 					EEG.event(index).type = 'topresp';
% $$$ 				else
% $$$ 					EEG.event(index).type = 'bottomresp';
% $$$ 				end;
% $$$ 				EEG.event(end+1).latency = EEG.event(index).latency;
% $$$ 				if (resultcode == targetcode) 
% $$$ 					EEG.event(end).type = 'correct';
% $$$ 				else
% $$$ 					EEG.event(end).type = 'miss';
% $$$ 				end;
% $$$ 			end;
% $$$ 		end;
% $$$ 	end;
% $$$ 	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
end;

command = sprintf('EEG = pop_loadbci(''%s'', %f);',filename, srate); 
return;
