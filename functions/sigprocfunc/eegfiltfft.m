% eegfiltfft() -  (high|low|band)-pass filter data using inverse fft 
%                 (without using the Matlab signal processing toolbox)
% Usage:
%  >> [smoothdata] = eegfiltfft(data,srate,locutoff,hicutoff);
%  >> [smoothdata] = eegfiltfft(data,srate,locutoff,hicutoff,epochframes,filtorder,revfilt);
%
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {0 -> lowpass}
%   hicutoff    = high-edge frequency in pass band (Hz) {0 -> highpass}
%
% Optional inputs:
%   epochframes = frames per epoch (filter each epoch separately {def/0: data is 1 epoch}
%   filtorder   = argument not used (but required for symetry with eegfilt() function).
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch filter). {0}
%
% Outputs:
%    smoothdata = smoothed data
%
% Known problems:
%    The signal drop off is much smaller compared to standard filtering methods
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2003
%
% See also: eegfilt()

% inspired from a ggogle group message
% http://groups.google.com/groups?q=without+%22the+signal+processing+toolbox%22&hl=en&lr=&ie=UTF-8&oe=UTF-8&selm=f56893ae.0311141025.3069d4f8%40posting.google.com&rnum=8

% Copyright (C) Arnaud Delorme, SCCN/INC/UCSD, arno@salk.edu
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

function smoothdata = eegfiltfft(data, fs, lowcut, highcut, epochframes, filtorder, revfilt);
    
    if nargin < 4
        help eegfiltfft;
    end;
    [chans frames] = size(data);
    if nargin < 5 | epochframes == 0
        epochframes = frames;
    end
    if nargin < 7
        revfilt = 0;
    end;
    
    epochs = frames/epochframes;
    if epochs ~= round(epochs)
        error('Epochframes does not divide the total number of frames');
    end;
    fv=reshape([0:epochframes-1]*fs/epochframes,epochframes,1); % Frequency vector for plotting    
    
    %figure;
    %plot(fv, 20*log(abs(X))/log(10))  % Plot power spectrum in dB scale
    %xlabel('Frequency [Hz]')
    %ylabel('Signal power [dB]')
    
    % find closest freq in fft decomposition
    % --------------------------------------
    if lowcut ~= 0
        [tmp idxl]=min(abs(fv-lowcut));  % Find the entry in fv closest to 5 kHz
    else
        idxl = 0;
    end;
    if highcut ~= 0        
        [tmp idxh]=min(abs(fv-highcut));  % Find the entry in fv closest to 5 kHz    
    else 
        idxh = ceil(length(fv)/2);
    end;
    
    % filter the data
    % ---------------
    smoothdata = zeros(chans,frames);
    for e = 1:epochs                % filter each epoch, channel 
        for c=1:chans
            X=fft(data(c,(e-1)*epochframes+1:e*epochframes));
            if revfilt
                X(idxl+1:idxh-1)=0;
                if mod(length(X),2) == 0
                    X(end/2:end)=0;
                else
                    X((end+1)/2:end)=0;
                end;
            else
                X(1:idxl)=complex(0);
                X(end-idxl:end)=complex(0);
                X(idxh:end)=complex(0);
            end;                
            smoothdata(c,(e-1)*epochframes+1:e*epochframes) = 2*real(ifft(X));
            if epochs == 1 
                if rem(c,20) ~= 0
                    fprintf('.');
                else 
                    fprintf('%d',c);
                end
            end
        end
    end
    fprintf('\n');

