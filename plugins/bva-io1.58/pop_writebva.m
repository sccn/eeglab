% pop_writebva() - export EEG dataset
% 
% Usage:
%   >> EEG = pop_writebva(EEG);   % a window pops up
%   >> EEG = pop_writebva(EEG, filename);
%
% Inputs:
%   EEG            - eeglab dataset
%   filename       - file name
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2005-

% Copyright (C) 2005, Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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

function com = pop_writebva(EEG, filename); 

com = '';
if nargin < 1 
    help pop_writebva;
    return;
end;

if nargin < 2
    [filename, filepath] = uiputfile('*', 'Output file');
    if length( filepath ) == 0 return; end;
    filename = [ filepath filename ];
end;

% remove extension if any
% -----------------------
posdot = find(filename == '.');
if ~isempty(posdot), filename = filename(1:posdot(end)-1); end;

% open output file
% ----------------
fid1 = fopen( [ filename '.vhdr' ], 'w' );
fid2 = fopen( [ filename '.vmrk' ], 'w' );
fid3 = fopen( [ filename '.dat'  ], 'wb', 'ieee-le');
[ tmppath basename ] = fileparts( filename );

% write data
% ----------
for index = 1:EEG.nbchan
    fwrite(fid3, EEG.data(index,:), 'float' );
end;

% write header
% ------------
fprintf(fid1, 'Brain Vision Data Exchange Header File Version 1.0\n');
fprintf(fid1, '; Data created from the EEGLAB software\n');
fprintf(fid1, '\n');
fprintf(fid1, '[Common Infos]\n');
fprintf(fid1, 'DataFile=%s\n', [ basename '.dat'  ]);
if ~isempty(EEG.event)
    fprintf(fid1, 'MarkerFile=%s\n', [ basename '.vmrk' ]);
end;
fprintf(fid1, 'DataFormat=BINARY\n');
fprintf(fid1, '; Data orientation: VECTORIZED=ch1,pt1, ch1,pt2..., MULTIPLEXED=ch1,pt1, ch2,pt1 ...\n');
fprintf(fid1, 'DataOrientation=VECTORIZED\n');
fprintf(fid1, 'DataType=TIMEDOMAIN\n');
fprintf(fid1, 'NumberOfChannels=%d\n', EEG.nbchan);
fprintf(fid1, 'DataPoints=%d\n', EEG.pnts*EEG.trials);
fprintf(fid1, '; Sampling interval in microseconds if time domain (convert to Hertz:\n');
fprintf(fid1, '; 1000000 / SamplingInterval) or in Hertz if frequency domain:\n');
fprintf(fid1, 'SamplingInterval=%d\n', 1000000/EEG.srate);
if EEG.trials > 1
    fprintf(fid1, 'SegmentationType=MARKERBASED\n');
end;
fprintf(fid1, '\n');
fprintf(fid1, '[Binary Infos]\n');
fprintf(fid1, 'BinaryFormat=IEEE_FLOAT_32\n');
fprintf(fid1, '\n');
if ~isempty(EEG.chanlocs)
    fprintf(fid1, '[Channel Infos]\n');
    fprintf(fid1, '; Each entry: Ch<Channel number>=<Name>,<Reference channel name>,\n');
    fprintf(fid1, '; <Resolution in microvolts>,<Future extensions..\n');
    fprintf(fid1, '; Fields are delimited by commas, some fields might be omited (empty).\n');
    fprintf(fid1, '; Commas in channel names are coded as "\1".\n');
    for index = 1:EEG.nbchan
        fprintf(fid1, 'Ch%d=%s,, \n', index, EEG.chanlocs(index).labels);
    end;
    fprintf(fid1, '\n');

    disp('Warning: channel location were not exported to BVA (it will use default');
    disp('         10-20 BESA locations based on channel names)');
    %if isfield(EEG.chanlocs, 'sph_radius')
    %    fprintf(fid1, '[Coordinates]\n');
    %    fprintf(fid1, '; Each entry: Ch<Channel number>=<Radius>,<Theta>,<Phi>\n');
    %    loc = convertlocs(EEG.chanlocs, 'sph2sphbesa');
    %    for index = 1:EEG.nbchan
    %        fprintf(fid1, 'Ch%d=%d,%d,%d\n', index, round(loc(index).sph_theta_besa), ...
    %                                        round(loc(index).sph_phi_besa), 0);
    %    end;
    %end;
end;

% export event information
% ------------------------
if ~isempty(EEG.event)
    fprintf(fid2, 'Brain Vision Data Exchange Marker File, Version 1.0\n');
    fprintf(fid2, '; Data created from the EEGLAB software\n');
    fprintf(fid2, '; The channel numbers are related to the channels in the exported file.\n');
    fprintf(fid2, '\n');
    fprintf(fid2, '[Common Infos]\n');
    fprintf(fid2, 'DataFile=%s\n', [ basename '.dat'  ]);
    fprintf(fid2, '\n');
    fprintf(fid2, '[Marker Infos]\n');
    fprintf(fid2, '; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,\n');
    fprintf(fid2, '; <Size in data points>, <Channel number (0 = marker is related to all channels)>,\n');
    fprintf(fid2, '; <Date (YYYYMMDDhhmmssuuuuuu)>\n');
    fprintf(fid2, '; Fields are delimited by commas, some fields might be omited (empty).\n');
    fprintf(fid2, '; Commas in type or description text are coded as "\1".\n');
    
    % rename type and comments
    % ------------------------
    for index = 1:length(EEG.event)
        EEG.event(index).comment = EEG.event(index).type;
        EEG.event(index).type    = 'Stimulus';
    end;
    
    % make event cell array
    % ---------------------
    for index = 1:EEG.trials
        EEG.event(end+1).type    = 'New Segment';
        EEG.event(end  ).latency = (index-1)*EEG.pnts+1;
    end;
    tmpevent = EEG.event;
    [tmp latorder ] = sort( [ tmpevent.latency ] );
    EEG.event = EEG.event(latorder);

    % Recode boundary events
    % ----------------------
    bndArray = find(strcmp('boundary', {tmpevent.comment})); % Find boundary events
    isDupArray = ismember([tmpevent(bndArray).latency], [tmpevent(strcmp('New Segment', {tmpevent.type})).latency]); % Find already existing New Segment events with identical latency
    notduplist  = bndArray(~isDupArray);
    duplist     = bndArray(isDupArray);
    for index = 1:length(notduplist)
        EEG.event(notduplist(index)).type    = 'New Segment'; % Recode boundary event type
        EEG.event(notduplist(index)).comment = '';            % Recode boundary event comment
    end;
    EEG.event(duplist) = []; % Remove duplicate New Segment events 

    % rename latency events
    % ---------------------
    time0ind = [];
    for index = 1:length(EEG.event)
        if mod( EEG.event(index).latency, EEG.pnts) == -EEG.xmin*EEG.srate+1
            time0ind = [ time0ind index ];
        end;
    end;
    for index = length(time0ind):-1:1
        EEG.event(time0ind(index)+1:end+1) = EEG.event(time0ind(index):end);
        EEG.event(time0ind(index)).type    = 'Time 0';
        EEG.event(time0ind(index)).comment = ''; 
    end;
    
    % write events
    % ------------
    e = EEG.event;
    for index = 1:length(e)
        
        % duration field
        % --------------
        if isfield(e, 'duration')
            if ~isempty(e(index).duration)
                tmpdur = e(index).duration;
            else tmpdur = 0;
            end;
        else tmpdur = 0;
        end;

        % comment field
        % -------------
        if isfield(e, 'comment')
            if ~isempty(e(index).comment)
                tmpcom = e(index).comment;
            else tmpcom = '';
            end;
        else tmpcom = num2str(e(index).type);
        end;

        fprintf(fid2, 'Mk%d=%s,%s,%d,%d,0,0\n', index, num2str(e(index).type), num2str(tmpcom), round(e(index).latency), tmpdur);
    end;
end
fclose(fid1);
fclose(fid2);
fclose(fid3);

com = sprintf('pop_writebva(%s,''%s'');', inputname(1), filename); 
return;
