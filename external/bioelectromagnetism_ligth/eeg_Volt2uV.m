function [uV] = eeg_volt2uv(Volts);

% eeg_volt2uv - Convert volts to microvolts
%
%   [uV] = eeg_volts2uv(Volts)
%
%   Simply, uV = Volts .* 10^6
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:50 $

% Licence:  GNU GPL, no implied or express warranties
% History:  11/2001, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uV = Volts .* 10^6;
