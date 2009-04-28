function coh = eeg_load_scan_coh_asc(filename)

% eeg_load_scan_coh_asc - Read a Scan coherence ascii file (*.dat)
% 
% Usage: coh = eeg_load_scan_coh_asc(filename)
% 
% This script extracts all the coherence values from a Scan 
% coherence ascii file.
% 
% An example of the coh struct returned:
% 
%       subject: ''
%          date: '19-Jul-2002'
%          time: '14:01:47'
%      channels: 64
%          rate: 500
%          type: 'PairedChannels'
%        points: 82
%        domain: 'Frequency'
%         pairs: 2016
%      freqbins: [1x82 double]
%     pairnames: {2016x2 cell}
%          data: [2016x82 double]
%         phase: [2016x82 double]
% 
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no express or implied warranties
% History:  11/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[path,name,ext] = fileparts(filename);
file = fullfile(path,[name ext]);

[fid,msg] = fopen(file,'r');
if ~isempty(msg), error(msg); end

tic

fprintf('\nEEG_LOAD_SCAN_COH_ASC...\n');
fprintf('...reading data file:\n\t%s\n',file);

coh = read_coh(fid);

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coh] = read_coh(fid)
    
    coh = [];
    
    % Read the header, eg:
	% 
	% [Subject]	
	% [Date]	07/19/2002
	% [Time]	14:01:47
	% [Channels]	64
	% [Rate]	500.000000
	% [Type]	PairedChannels
	% [Points]	82
	% [Domain]	Frequency
	% [Pairs]	2016
	% [Paired Channels Data Matrix: Magnitude Squared (Pairs by Points)]
    
    fprintf('...reading header\n');
    
    while 1,
        tmp = fgetl(fid);
        if strmatch(lower('[Subject]'),lower(tmp)),
            coh.subject = strrep(tmp,sprintf('[Subject]\t'),'');
            tmp = fgetl(fid);
        end
        if strmatch(lower('[Date]'),lower(tmp)),
            tmp = strrep(tmp,sprintf('[Date]\t'),'');
            coh.date = datestr(datenum(tmp));
            tmp = fgetl(fid);
        end
        if strmatch(lower('[Time]'),lower(tmp)),
            coh.time = strrep(tmp,sprintf('[Time]\t'),'');
            tmp = fgetl(fid);
        end
        if strmatch(lower('[Channels]'),lower(tmp)),
            tmp = strrep(tmp,sprintf('[Channels]\t'),'');
            coh.channels = str2num(tmp);
            tmp = fgetl(fid);
        end
        if strmatch(lower('[Rate]'),lower(tmp)),
            tmp = strrep(tmp,sprintf('[Rate]\t'),'');
            coh.rate = str2num(tmp);
            tmp = fgetl(fid);
        end
        if strmatch(lower('[Type]'),lower(tmp)),
            coh.type = strrep(tmp,sprintf('[Type]\t'),'');
            tmp = fgetl(fid);
        end
        if strmatch(lower('[Points]'),lower(tmp)),
            tmp = strrep(tmp,sprintf('[Points]\t'),'');
            coh.points = str2num(tmp);
            tmp = fgetl(fid);
        end
        if strmatch(lower('[Domain]'),lower(tmp)),
            coh.domain = strrep(tmp,sprintf('[Domain]\t'),'');
            tmp = fgetl(fid);
        end
        if strmatch(lower('[Pairs]'),lower(tmp)),
            tmp = strrep(tmp,sprintf('[Pairs]\t'),'');
            coh.pairs = str2num(tmp);
            tmp = fgetl(fid);
        end
        if findstr(lower(tmp),lower('Data')),
            datamatrix = tmp;
            tmp = fgetl(fid);
            tmp = strrep(tmp,'[' ,'');
            tmp = strrep(tmp,']' ,'');
            tmp = strrep(tmp,'Hz','');
            coh.freqbins = str2num(tmp);
        end
        break;
    end
    
    coh.pairnames = cell(coh.pairs,2);
    
    if findstr(datamatrix,'Pairs by Points'),
        coh.data = zeros(coh.pairs,coh.points);
        row = coh.pairs;
    elseif findstr(datamatrix,'Points by Pairs'),
        % not sure this is an option, but just in case
        fprintf('..cannot read data matrix for points by pairs\n');
        return
        %coh.data = zeros(coh.points,coh.pairs);
        %row = coh.points;
    end
    
    % Now read in the data matrix
    fprintf('...reading coherence\n');
    for r = 1:row,
        tmp = fgetl(fid);
        a = findstr(tmp,'[') + 1;
        b = findstr(tmp,']') - 1;
        c = findstr(tmp,'-');
        c = c(1);
        coh.pairnames{r,1} = tmp(a:c-1);
        coh.pairnames{r,2} = tmp(c+1:b);
        coh.data(r,:) = str2num(tmp(b+2:end));
    end
    
    % Phase Degrees
    tmp = fgetl(fid);
    if findstr(tmp,'Phase'),
        fprintf('...reading phase (degrees)\n');
        tmp = fgetl(fid); % this line contains freqbins
        coh.phase = zeros(size(coh.data));
        for r = 1:row,
            tmp = fgetl(fid);
            b = findstr(tmp,']') - 1;
            coh.phase(r,:) = str2num(tmp(b+2:end));
        end
    end
    
    fclose(fid);
    
return
