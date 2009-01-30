function ctf_print_setup(ctf);

% ctf_print_setup - print setup information to matlab command window
%
% ctf_print_setup(ctf)
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                      >
%      <                    DISCLAIMER:                       >
%      <                                                      >
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                    OFFICIAL USE.                     >
%      <                                                      >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
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

% Modified: 02/2004, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('...data format       : %s\n',ctf.res4.header);
fprintf('...date collected    : %s\n',ctf.setup.date);
fprintf('...time collected    : %s\n',ctf.setup.time);
fprintf('...run name          : %s\n',ctf.setup.run_name);
fprintf('...run title         : %s\n',ctf.setup.run_title);
fprintf('...subject           : %s\n',ctf.setup.subject);
fprintf('...run description   : %s\n',ctf.setup.run_description);
fprintf('...operator          : %s\n',ctf.setup.operator);
fprintf('...number of channels: %d\n',ctf.setup.number_channels);
fprintf('...number of samples : %d per trial\n',ctf.setup.number_samples);
fprintf('...sample rate       : %g samples/sec\n',ctf.setup.sample_rate);
fprintf('...number of trials  : %g trials (average of %g)\n',ctf.setup.number_trials,ctf.setup.number_trials_averaged);
fprintf('...duration          : %g sec total\n',ctf.setup.duration_total);
fprintf('...duration          : %g sec / trial\n',ctf.setup.duration_trial);
fprintf('...pretrigger points : %g samples\n',ctf.setup.pretrigger_samples);
fprintf('...pretrigger msec   : %g msec\n',ctf.setup.pretrigger_msec);
fprintf('...sensor file name  : %s\n',ctf.setup.sensor_file_name);
fprintf('...head zeroing      : %s\n',ctf.setup.head_zero);
fprintf('...number of filters : %d\n',ctf.setup.number_filters);
for i = 1:ctf.setup.number_filters,
  fprintf(' -Filter # \t%g\n',i)
  fprintf(' -Frequency: \t%g Hz\n',ctf.setup.filters(i).freq)
  fprintf(' -Class: \t%g\n',ctf.setup.filters(i).class)
  fprintf(' -Type: \t%g\n',ctf.setup.filters(i).type)
  if ~isempty(ctf.setup.filters(i).params)
    fprintf(' -Parameter(s): \t%g\n',ctf.setup.filters(i).params)
  end
end
%fprintf('\n');

return


% from ctf_read_res4

% for i = 1:ctf.setup.number_filters,
%   ctf.setup.filters(i).freq     = fread(fid,1,'double');
%   ctf.setup.filters(i).class    = fread(fid,1,'int32');
%   ctf.setup.filters(i).type     = fread(fid,1,'int32');
%   ctf.setup.filters(i).numparam = fread(fid,1,'int16');
%   ctf.setup.filters(i).params   = fread(fid,ctf.setup.filters(i).numparam,'double');
% end
