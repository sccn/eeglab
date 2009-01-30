function ctf = ctf_filter(ctf,lcf,hcf,order);

% ctf_filter - apply a butterworth polynomial filter
% 
% 	Usage	: ctf = ctf_filter(ctf,lcf,hcf,order);
%
% 	- input arguments 
%		ctf   : meg data file returned by ctf_read
%		lcf   : low cutoff frequency (default 0.01)
%		hcf	  : high cutoff frequency (default 40)
%		order : butterworth polynomial order (default 2)
%
% 	- output argument
%		ctf.data : filtered replacement of ctf.data
%
%  This function calls the butter and filtfilt functions of the matlab
%  signal processing toolbox.  The filtfilt values are baselined with
%  ctf_baseline, so the baseline offsets can be slightly different from the
%  input ctf.data.
%


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2004  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Created: 05/2004, copyright 2004 Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% low cutoff frequency (default 2)
if ~exist('lcf','var') || isempty(lcf),
    lcf = 0.01;
end
% high cutoff frequency (default 20)
if ~exist('hcf','var') || isempty(hcf),
    hcf = 40;
end
% butter filter order (default 2)
if ~exist('order','var') || isempty(order),
    order = 2;
end

sample_freq = ctf.setup.sample_rate;

if hcf > (sample_freq/2),
    warning('hcf > sample_freq/2, setting hcf = sample_freq/2');
    hcf = sample_freq / 2;
end
if lcf <= 0 || lcf > (sample_freq/2) || lcf >= hcf,
    warning('lcf value is <=0 or >(sample_freq/2) or >=hcf, setting lcf = 2');
    lcf = 2;
end


% design the Butterworth filter
cf1 = lcf/(sample_freq/2);
cf2 = hcf/(sample_freq/2);
[B,A] = butter(order,[cf1 cf2]);

% filter the data
% data should be N samples x M channels
for trial = 1:size(ctf.data,3)
    data = ctf.data(:,:,trial);
    ctf.data(:,:,trial) = filtfilt(B,A,data);
end

% now rebaseline the filtered data
ctf = ctf_baseline(ctf);

return
