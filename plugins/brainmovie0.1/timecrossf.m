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
%        = timecrossf(data, frames, tlimits, srate, cycles, ...
%                 winsize,timesout, padratio,maxfreq, tvec,eloc_file, ...
%                 alpha,marktimes,powbase,pboot,rboot);
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

% $Log: not supported by cvs2svn $
% Revision 1.2  2002/11/18 00:53:43  arno
% new version
%
% Revision 1.1  2002/11/18 00:43:32  arno
% Initial revision
%
% Revision 1.1  2002/11/18 00:33:30  arno
% Initial revision
%

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
end;

% compute all time frequency maps
% -------------------------------
compter = 1;
for numcompo = 1:nbcompo
	if iscell(data)
        [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef({ data1(numcompo,:) data2(numcompo,:) }, ...
                                                      frames, tlimits, srate, cycle, varargin{:});   
        ALLERSP1{numcompo,1}     = applyboot( ersp{1}, erspboot{1});
        ALLERSP2{numcompo,1}     = applyboot( ersp{2}, erspboot{2});
        ALLERSPDIFF{numcompo,1}  = applyboot( ersp{3}, erspboot{3});
        ALLITC1{numcompo,1}      = applyboot( abs(itc{1}) , itcboot{1});
        ALLITC2{numcompo,1}      = applyboot( abs(itc{2}) , itcboot{2});
        ALLITCDIFF{numcompo,1}   = applyboot( abs(itc{3}) , itcboot{3});
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
                ALLCROSSF1     { index1, index2 } = applyboot(abs(coh{1}), cohboot{1});
                ALLCROSSF2     { index1, index2 } = applyboot(abs(coh{2}), cohboot{2});
                ALLCROSSFDIFF  { index1, index2 } = applyboot(abs(coh{3}), cohboot{3});
                ALLCROSSFANGLE1{ index1, index2 } = cohangles{1};
                ALLCROSSFANGLE2{ index1, index2 } = cohangles{2};
                ALLCROSSFANGLEDIFF{ index1, index2 } = cohangles{3};
            else
                [coh,mcoh,timesout,freqsout,cohboot,cohangles] = newcrossf(data(index1,:), frames,  ...
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
if iscell(data)
    ALLERSP    = { ALLERSP1    ALLERSP2    ALLERSPDIFF };
    ALLITC     = { ALLITC1     ALLITC2     ALLITCDIFF  };
    ALLCROSSF  = { ALLCROSSF1  ALLCROSSF2  ALLCROSSFDIFF  };
    ALLCROSSFANGLE  = { ALLCROSSFANGLE1  ALLCROSSFANGLE2  ALLCROSSFANGLEDIFF  };
end;

return;

function array = applyboot(array, arrayboot)
    if size(erspboot,3) == 2
        array(find((array > arrayboot(:,:,1)) & (array < arrayboot(:,:,2)))) = 0;
    elseif size(erspboot,2) ~= 2
        array(find(array > arrayboot(:,:))) = 0;
    else        
        array(find((array > repmat(arrayboot(:,1),[1 size(array,2)])) & ...
                   (array < repmat(arrayboot(:,2),[1 size(ersp,2)])))) = 0;
    end;
