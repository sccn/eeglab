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

function [EEG, command] = pop_loadbci(filename, srate); 

    EEG = [];
    command = '';
    
    if nargin < 1
        % ask user
        [filename, filepath] = uigetfile('*.*', 'Choose a BCI file -- pop_loadbci'); 
        drawnow;
        if filename == 0 return; end;
        filename = [filepath filename];
        promptstr    = { 'Sampling rate' };
        inistr       = { '256' };
        result       = inputdlg2( promptstr, 'Import BCI2000 data -- pop_loadbci()', 1,  inistr, 'pop_loadbci');
        if length(result) == 0 return; end;
        srate   = eval( result{1} );
    end;
    
    % import data
    % -----------
    EEG = eeg_emptyset;
    fprintf('Pop_loadbci: importing BCI file...\n');
    try
        % try to read as matlab
        % ---------------------
        bci = load( filename, '-mat');
        allfields = fieldnames(bci);
        allfields = setdiff(allfields, 'signal');
        for index = 1:size(bci.signal,2)
            chanlabels{index} = [ 'C' int2str(index) ];
        end;
        for index = 1:length(allfields)
            bci.signal(:,end+1) = getfield(bci, allfields{index});
        end;
        EEG.chanlocs = struct('labels', { chanlabels{:} allfields{:} });
        EEG.data     = bci.signal';
        EEG.nbchan   = size(EEG.data, 1);
        EEG.pnts     = size(EEG.data, 2);
        EEG.trials   = 1;
        EEG.srate    = srate;
        EEG.comments = [ 'Original file: ' filename ];
        EEG = eeg_checkset(EEG);
        return;
    catch
        % get file names
        % --------------
        fields = loadtxt(filename, 'nlines', 1, 'verbose', 'off');
        if length(fields) > 300
            error('Not a BCI ASCII file');
        end;
        
        % read data
        % ---------
        fid = fopen(filename, 'r');
        allcollumns = fgetl(fid);
        tmpdata = fscanf(fid, '%d', Inf);
        tmpdata = reshape(tmpdata, length(fields), length(tmpdata)/length(fields));
        
        EEG.data = tmpdata;
        EEG.chanlocs = struct('labels', fields);
        EEG.nbchan = size(EEG.data, 1);
        EEG.pnts   = size(EEG.data, 2);
        EEG.trials = 1;
        EEG.srate  = srate;
        EEG = eeg_checkset(EEG);
        return;

        % data channel range
        % ------------------
        indices = strmatch('ch', fields);
        
        bci = [];
        for index = setdiff(1:length(fields), indices)
            bci = setfield(bci, fields{index}, tmpdata(index,:));
        end;
        bci.signal = tmpdata(indices,:);
    end;
    
    
    % ask for which event to import
    % -----------------------------
    geom = {[0.7 0.7 0.7]};
    uilist = { { 'style' 'text' 'string' 'State name' 'fontweight' 'bold' } ...
               { 'style' 'text' 'string' '    Import' 'fontweight' 'bold'  } ...
               { 'style' 'text' 'string' 'Type of'  'fontweight' 'bold' } };
    allfields = setdiff(fieldnames(bci), 'signal');
    latencyfields = { '-----' };
    for index = 1:length(allfields)
        if ~isempty(findstr( lower(allfields{index}), 'time')) 
            latencyfields{end+1} = allfields{index};
        end;
    end;
    for index = 1:length(allfields)
        geom   = { geom{:} [1.3 0.3 0.3 0.3 1] };
        uilist{end+1} = { 'style' 'text' 'string' allfields{index} };
        if ~isempty(findstr( lower(allfields{index}), 'time'))
            uilist{end+1} = { 'style' 'checkbox' 'value' 0 };
            uilist{end+1} = { };
            uilist{end+1} = { };
            uilist{end+1} = { };
        else
            uilist{end+1} = { 'style' 'checkbox' };
            uilist{end+1} = { };
            uilist{end+1} = { };
            uilist{end+1} = { 'style' 'listbox' 'string' strvcat(latencyfields) };
        end;
    end;
    geom   = { geom{:} [1] [0.08 1] };
    uilist{end+1} = { };
    uilist{end+1} = { 'style' 'checkbox' 'value' 0 };
    uilist{end+1} = { 'style' 'text' 'string' 'Attempt to adjust event latencies using sourcetime?' };
    
    result = inputgui( geom, uilist, 'pophelp(''pop_loadbci'')', 'Import BCI2000 data files - pop_loadbci()');
    if isempty(result), return; end;
    
    % convert results to command line input
    % -------------------------------------
    listimport = {};
    count = 1;
    for index = 1:length(allfields)
        if ~isempty(findstr( lower(allfields{index}), 'time')) 
            if result{count}, listimport{end+1} = 'event'; listimport{end+1} = { allfields{index} }; end;
            count = count+1;
        else 
            if result{count}
                if result{count+1} ~= 1
                    listimport{end+1} = 'event'; listimport{end+1} = { allfields{index}  allfields{result{count+1}-1} };
                else
                    listimport{end+1} = 'event'; listimport{end+1} = { allfields{index} };
                end;
            end;
            count = count+2;                
        end;
    end;
    if result{end}, adjust = 1; else adjust = 0; end;
    
    % decode command line input 
    % -------------------------
    count = 1;
    for index = 2:2:length(listimport)
        tmpindmatch = strmatch(listimport{index}{1}, allfields, 'exact');
        if ~isempty(tmpindmatch), indeximport(count) = tmpindmatch; 
        else error(['State ''' listimport{index}{1} ''' not found']); 
        end;
        if length( listimport{index} ) > 1
            tmpindmatch = strmatch(listimport{index}{2}, allfields, 'exact');
            if ~isempty(tmpindmatch), corresp(count) = tmpindmatch; 
            else error(['State ''' listimport{index}{2} ''' not found']); 
            end;
        else
            corresp(count) = 0;
        end;
        count = count+1;
    end;
        
    % find block size
    % ---------------
    tmpevent = find( diff(getfield(bci, 'SourceTime')) ~= 0);
    diffevent = tmpevent(2:end)-tmpevent(1:end-1);
    blocksize = unique(diffevent);
    if length(blocksize) > 1, error('Error in determining block size'); 
    else                      fprintf('Blocksize: %d\n', blocksize); 
    end;
    
    % find types
    % ----------
    tmpcorresp      = find(corresp);
    indexcorresp    = corresp(tmpcorresp);
    indexcorrespval = indeximport(tmpcorresp);
    if length(tmpcorresp) ~= length(intersect(corresp, indeximport))
        disp('Warning: correspondance problem, some information will be lost');
    end;
    
    % remove type from latency array
    % ------------------------------
    indeximport(tmpcorresp) = [];
    if adjust
        fprintf('Latency of event will be adjusted\n');
    else fprintf('WARNING: Latency of event will not be adjusted (latency uncertainty %2.1f ms)\n', ...
                 blocksize/srate*1000);
    end;
    
	% process events
	% --------------
	fprintf('Pop_loadbci: importing events...\n');
	counte = 1; % event counter
	events(10000).latency = 0;
    indexsource =  strmatch('sourcetime', lower( allfields ), 'exact' );
    sourcetime  = getfield(bci, allfields{ indexsource });
	for index = 1:length(indeximport)
        tmpdata  = getfield(bci, allfields{indeximport(index)});
		tmpevent = find( diff(tmpdata) > 0);
		tmpevent = tmpevent+1;
        tmpcorresp = find(indexcorresp == indeximport(index));
        
		for tmpi = tmpevent'
            if ~isempty(tmpcorresp)
                events(counte).type = allfields{ indexcorrespval(tmpcorresp) };
            else 
                events(counte).type = allfields{ indeximport(index) };
            end;
            if adjust
                baselatency = sourcetime(tmpi); % note that this is the first bin a block
                realtmpi    = tmpi+blocksize;   % jump to the end of the block+1
                if curlatency < baselatency, curlatency = curlatency+65536; end; % in ms
                events(counte).latency = realtmpi+(curlatency-baselatency)/1000*srate;
                % (curlatency-baselatency)/1000*srate
                % there is still a potentially large error between baselatency <-> realtmpi
            else
                events(counte).latency = tmpi+(blocksize-1)/2;
            end;
            counte = counte+1;
		end;
	end; 
           
	EEG.data = bci.signal';
	EEG.nbchan = size(EEG.data, 1);
	EEG.pnts   = size(EEG.data, 2);
	EEG.trials = 1;
	EEG.srate  = srate;
	EEG = eeg_checkset(EEG);
	EEG.event = events(1:counte-1);	
	EEG = eeg_checkset( EEG, 'eventconsistency' );
	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
    
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ dsffd
% $$$         
% $$$ 	% find electrode indices
% $$$ 	% ----------------------
% $$$ 	indices = [];
% $$$ 	for index = 1:length(colnames)
% $$$ 		if strcmp( colnames{index}(1:2), 'ch')
% $$$ 			indices = [indices index ];
% $$$ 		end;
% $$$ 	end;
% $$$     
% $$$     EEG.data = tmpdata(indices,:);
% $$$     EEG.nbchan = size(EEG.data, 1);
% $$$     EEG.srate  = srate;
% $$$     try
% $$$         eventindices = setdiff(1:length(colnames), indices);
% $$$         ISIind = eventindices(3 + 9);
% $$$         eventindices(3 + [ 1 2 3 4 7 8 9 10 11 12]) = [];
% $$$         eventindices(1:3) = []; % suppress these event 
% $$$         
% $$$         % add the trial number
% $$$         % --------------------
% $$$         tmptrial = find( diff(tmpdata(ISIind, :)) ~= 0);
% $$$         tmptrial = tmptrial+1;
% $$$         
% $$$         % process events
% $$$         % --------------
% $$$         fprintf('Pop_loadbci: importing events...\n');
% $$$         counte = 1; % event counter
% $$$         events(10000).latency = 0;
% $$$         for index = eventindices
% $$$             counttrial = 1;
% $$$             tmpevent = find( diff(tmpdata(index, :)) ~= 0);
% $$$             tmpevent = tmpevent+1;
% $$$             for tmpi = tmpevent
% $$$                 if tmpdata(index, tmpi)
% $$$                     events(counte).type    = [ colnames{index} int2str(tmpdata(index, tmpi)) ];
% $$$                     events(counte).latency = tmpi;
% $$$                     %events(counte).value   = tmpdata(index, tmpi);
% $$$                     %while tmpi > tmptrial(counttrial) & counttrial < length(tmptrial)
% $$$                     %	counttrial = counttrial+1;
% $$$                     %end;
% $$$                     %events(counte).trial = counttrial;				
% $$$                     counte = counte+1;
% $$$                     %if mod(counte, 100) == 0, fprintf('%d ', counte); end;
% $$$                 end;
% $$$             end;
% $$$         end;
% $$$ 	
% $$$         % add up or down events
% $$$         % ---------------------
% $$$         EEG = eeg_checkset(EEG);
% $$$         EEG.event = events(1:counte-1);	
% $$$         EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
% $$$         for index=1:length(EEG.event)
% $$$             if strcmp(EEG.event(index).type(1:6), 'Target')
% $$$                 targetcode = str2num(EEG.event(index).type(end));
% $$$                 if targetcode == 1
% $$$                     EEG.event(index).type = 'toptarget';
% $$$                 else
% $$$                     EEG.event(index).type = 'bottomtarget';
% $$$                 end;
% $$$             else 
% $$$                 if strcmp(EEG.event(index).type(1:6), 'Result')
% $$$                     resultcode = str2num(EEG.event(index).type(end));
% $$$                     if resultcode == 1
% $$$                         EEG.event(index).type = 'topresp';
% $$$                     else
% $$$                         EEG.event(index).type = 'bottomresp';
% $$$                     end;
% $$$                     EEG.event(end+1).latency = EEG.event(index).latency;
% $$$                     if (resultcode == targetcode) 
% $$$                         EEG.event(end).type = 'correct';
% $$$                     else
% $$$                         EEG.event(end).type = 'miss';
% $$$                     end;
% $$$                 end;
% $$$             end;
% $$$         end;
% $$$         EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
% $$$         %EEG.data = tmpdata([72 73 75],:);
% $$$     catch, disp('Failed to import data events');
% $$$     end;
% $$$ 
% $$$ command = sprintf('EEG = pop_loadbci(''%s'', %f);',filename, srate); 
% $$$ return;
