function [CHAN1020,XYZ1020] = elec_1020select(CHAN)

% elec_1020select - select 10-20 locations
% 
% [labels,xyz] = elec_1020select(CHAN)
%
% where CHAN input is a cell array of channel names from the International
% 10-20 nomenclature for EEG electrode placement.  For a list of the 10-20
% electrode names, see the elec_1020all_cart function, which is based on:
%
% Oostenveld, R. & Praamstra, P. (2001). The five percent electrode system
% for high-resolution EEG and ERP measurements. Clinical Neurophysiology,
% 112:713-719.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Copyright (C) 2005  Darren L. Weber
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

% Modified: 01/2005, Darren.Weber_at_radiology.ucsf.edu
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $';
fprintf('\nELEC_1020SELECT [v %s]\n',ver(11:15));

% get the 1020 data
elec = elec_1020all_cart;
elec = struct2cell(elec);
labels = squeeze(elec(1,:,:))';
x = squeeze(elec(2,:,:)); x = x{1};
y = squeeze(elec(3,:,:)); y = y{1};
z = squeeze(elec(4,:,:)); z = z{1};
clear elec

% find all the electrode names in elec.labels that match CHAN
CHAN1020 = zeros(1,length(CHAN));
XYZ1020  = zeros(length(CHAN),3);
for c = 1:length(CHAN),
    chan = CHAN{c};
    index = strmatch(lower(chan),lower(labels),'exact');
    if ~isempty(index),
        CHAN1020(c) = index;
        XYZ1020(c,:) = [ x(index), y(index), z(index) ];
    else
        msg = sprintf('No match for channel: %s\n',chan)
        error(msg)
    end
end

CHAN1020 = labels(CHAN1020);

return
