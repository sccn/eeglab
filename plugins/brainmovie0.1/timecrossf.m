% timecrossf() - compute the time frequency content of an array rows and 
%                cross-coherence between the rows for. If the rows depict 
%                different components activation, these results can then 
%                be used to analyse time relation between these
%                components  (using the brainmovie function). You may 
%                also input electrode activities into this function to 
%                analyse the synchronisation and phase relashionship 
%                between the electrodes.  
%
% Usage:
%    >> [ allersps, allitcs, allcrossfs, allcrossfphis, times, freqs] ...
%        = timecrossf(data, frames, tlimits, srate, cycles, ...);
%
% Inputs:
%   data        - multi-row data (nb_input,frames*nepochs) for the first
%                 condition. Each row is a different input. The function 
%                 will conpute the time frequency content of each row and  
%                 the cross-coherence between the rows. If data is a cell
%                 array { data1 data2 }, the function will compute the 
%                 time-frequency decomposition for these 2 arrays and the 
%                 difference between them.
%   frames      - frames per epoch                        {768}
%   tlimits     - epoch time limits (ms) [mintime maxtime]{[-1000 2000]}
%   srate       - data sampling rate (Hz)                 {256}
%   ...         - see timef or crossf command for details about other 
%                 parameters 
%
% Outputs:
%   allersps      - cell array (size nb_row,1) of ERSP (event-related 
%                   power spectrum content of each row of the data array)
%   allitcs       - cell array (size nb_row,1) of ITC (event-related  
%                   inter-trial coherence content of each row of the data
%                   array)
%   allcrossfs    - cell array (size nb_row, nb_row) of cross-coherence 
%                   between the row of the data array (only the upper 
%                   diagonal part of the cell array is used)
%   allcrossfphis - cell array (size nb_row, nb_row) of cross-coherence 
%                   phase between the row of the data array (only the 
%                   upper diagonal part of the cell array is used)
%   times         - time array returned by the crossf or timef functions
%   freqs         - frequency array returned by the crossf or timef 
%                   functions
%
% Example1: generate a movie, 1 condition, below 10Hz dynamics  
%
%   % ICA decomposition, activity in array actICA (31 x (176 x n)),  
%   % graph for 6 components, 176 points from -100 ms to 600 ms after
%   % stimulus onset, 250 Hz sampling rate 
%   [ allersps, allitcs, allcrossfs, allcrossfphis, times, freqs] ...
%        = timecrossf( actICA(1:6,:), 176, [-100 600], 250, 1, 32, 100);
%   brainmovie( ersps, itcs, crossfs_amp, crossfs_phase, times, [1:2] )    
%   %[1:2] indicates the frequency rows to take into account (freqs(1:2))
%
% Example2: generate a movie, 2 condition, below 10Hz dynamics  
%
%   % Same as above with actICA2 (31 x (176 x n2)) the activity of the 
%   % components for a different condition
%   [ allersps, allitcs, allcrossfs, allcrossfphis, times, freqs] ...
%        = timecrossf( actICA(1:6,:), 176, [-100 600], 250, 1, 32, 100);
%   [ allersps(:,2), allitcs(:,2), allcrossfs(:,:,2), ...
%	                          allcrossfphis(:,:,2)] ...
%        = timecrossf( actICA2(1:6,:), 176, [-100 600], 250, 1, 32, 100);
%   brainmovie( allersps, allitcs, allcrossfs, allcrossfphis, times, ...
%                [1:2] )     
%
% See also timef(), crossf(), brainmovie()

% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or
% modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function [ ALLERSP, ALLITC, ALLCROSSF, ALLCROSSFANGLE, times, freqs ] = timecrossf(data, frames, tlimits, srate, cycle, varargin);

if nargin < 3
	help timecrossf;
	return;
end;

SAVE    = 0;	% save graphs in jpeg format
%XDISP   = 'cole:0.0';
nbcompo = size( data, 1);
if iscell(data)
    data1 = data{1};
    data2 = data{2};
    nbcompo = size( data1, 1);    
end;

% compute all time frequency maps
% -------------------------------
compter = 1;
for numcompo = 1:nbcompo
	if iscell(data)
        [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef({ data1(numcompo,:) data2(numcompo,:) }, ...
                                                          frames, tlimits, srate, cycle, varargin{:});
        size(ersp{1})
        ALLERSP{numcompo,1}     = applyboot( ersp{1}, erspboot{1});
        ALLERSP{numcompo,2}     = applyboot( ersp{2}, erspboot{2});
        ALLERSP{numcompo,3}     = applyboot( ersp{3}, erspboot{3});
        ALLITC {numcompo,1}     = applyboot( abs(itc{1}) , itcboot{1});
        ALLITC {numcompo,2}     = applyboot( abs(itc{2}) , itcboot{2});
        if ~isreal(itc{3}), itc{3} = abs(itc{3}); end;
        ALLITC {numcompo,3}     = applyboot( itc{3}, itcboot{3});
    else
        [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef( data(numcompo,:), ...
                                                          frames, tlimits, srate, cycle, varargin{:});   
        ALLERSP{numcompo,1}    = applyboot( ersp, erspboot);
        ALLITC{numcompo,1}     = applyboot( itc,  itcboot);
    end;
	if SAVE
		%command = sprintf('print -djpeg timef%2.2d', numcompo); eval(command); close;
		command = sprintf('hgsave(''timef%2.2d'');', numcompo); eval(command); close;
	end;
	compter = compter + 1 ;
end;	

% compute all cross coherence maps
% --------------------------------
ALLCROSSF = cell(nbcompo, nbcompo);
ALLCROSSFANGLE = cell(nbcompo, nbcompo);
for index1 = 1:nbcompo
	for index2 = 1:nbcompo
		if index2 > index1
            if iscell(data)
                [coh,mcoh,timesout,freqsout,cohboot,cohangles] = newcrossf({ data1(index1,:) data2(index1,:)}, ...
                                                                  { data1(index2,:) data2(index2,:)}, frames,  ...
                                                                  tlimits, srate, cycle, varargin{:});    
                size(coh{1})
                ALLCROSSF      { index1, index2, 1 } = applyboot(coh{1}, cohboot{1});
                ALLCROSSF      { index1, index2, 2 } = applyboot(coh{2}, cohboot{2});
                ALLCROSSF      { index1, index2, 3 } = applyboot(coh{3}, cohboot{3});
                ALLCROSSFANGLE { index1, index2, 1 } = cohangles{1};
                ALLCROSSFANGLE { index1, index2, 2 } = cohangles{2};
                ALLCROSSFANGLE { index1, index2, 3 } = cohangles{3};
            else
                [coh,mcoh,timesout,freqsout,cohboot,cohangles] = newcrossf(data(index1,:), data(index2,:), frames,  ...
                                                                  tlimits, srate, cycle, varargin{:});    
                ALLCROSSF      { index1, index2 } = applyboot(coh, cohboot);
                ALLCROSSFANGLE { index1, index2 } = cohangles;
                if SAVE
                    %command = sprintf('print -djpeg crossf%2.2d-%2.2d%', index1, index2); eval(command); close;
                    command = sprintf('hgsave(''crossf%2.2d-%2.2d%'');', index1, index2); eval(command); close;
                end;
            end;
		end;
	end;
end;

return;

function array = applyboot(array, arrayboot)
    if isempty(arrayboot), return; end;
    if size(arrayboot,3) == 2
        array(find((array > arrayboot(:,:,1)) & (array < arrayboot(:,:,2)))) = 0;
    elseif size(arrayboot,2) > 2
        array(find(array < arrayboot(:,:))) = 0;
    elseif size(arrayboot,2) == 2        
        array(find((array > repmat(arrayboot(:,1),[1 size(array,2)])) & ...
                   (array < repmat(arrayboot(:,2),[1 size(array,2)])))) = 0;
    else 
        array(find(array < repmat(arrayboot(:,1),[1 size(array,2)]))) = 0;    
    end;
