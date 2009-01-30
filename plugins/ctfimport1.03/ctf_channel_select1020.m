function [CHAN,ctf] = ctf_channel_select1020(ctf,CHAN,mod,plot)

% ctf_channel_select1020 - select CTF channels closest to 10-20 locations
% 
% [CHAN,ctf] = ctf_channel_select1020(ctf,CHAN,mod,plot)
%
% where CHAN input is a cell array of channel names from the International
% 10-20 nomenclature for EEG electrode placement.  The default is a basic
% 19 channel layout.  For a full list of 10-20 electrode names, see the
% elec_1020all_cart function, which is based on:
%
% Oostenveld, R. & Praamstra, P. (2001). The five percent electrode system
% for high-resolution EEG and ERP measurements. Clinical Neurophysiology,
% 112:713-719.
%
% Of course, the MEG sensors are not placed on head locations that move
% with the subject, as an EEG electrode will, so this function is an
% idealization.  This function fits the EEG/MEG sensors to a unit sphere
% and then finds the channels nearest the 10-20 locations requested.  The
% 10-20 locations are stored in the text file "elec_1020all_cart.txt" in
% this ctf toolbox.
%
% CHAN output is a numeric array of channel indices
% ctf output is modified to select only the 1020 channels if the mod input
% variable is 1.
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


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $

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

% Modified: 02/2004, Darren.Weber_at_radiology.ucsf.edu
%           03/2005, Darren.Weber_at_radiology.ucsf.edu
%                    added mod option, removed type return value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $';
fprintf('\nCTF_CHANNEL_SELECT1020 [v %s]\n',ver(11:15));

if ~exist('CHAN','var'), CHAN = []; end
if isempty(CHAN),
    % select 19 locations of 10-20 channel set
    CHAN = {'Fp1','Fp2',...
        'F7','F3','Fz','F4','F8',...
        'T7','C3','Cz','C4','T8',...
        'P7','P3','Pz','P4','P8',...
        'O1','O2'};
end

if ~exist('mod','var'), mod = 0; end
if isempty(mod), mod = 0; end
if ~exist('plot','var'), plot = 0; end
if isempty(plot), plot = 0; end

% get the 1020 data
[CHAN1020,XYZ1020] = elec_1020select(CHAN);

% translate the CTF sensors into a unit sphere, 
% assuming the origin is at (0,0,0)
sens = ctf.sensor.location';
distance = sqrt( sens(:,1).^2 + sens(:,2).^2 + sens(:,3).^2 );
distance = repmat(distance,1,3);
sens_unit = sens ./ distance;
clear distance

% K = DSEARCHN(X,T,XI) returns the indices K of the closest points in X for
%     each point in XI. X is an m-by-n matrix representing m points in n-D
%     space. XI is a p-by-n matrix, representing p points in n-D space.

CHAN = dsearchn(sens_unit,XYZ1020);
%[sortedCHAN,i] = unique(CHAN);


if plot
    figure;  hold on
    % plot all meg sensors
    scatter3( sens_unit(:,1),    sens_unit(:,2),    sens_unit(:,3), 'k');
    % plot the 10-20 locations required
    H1 = scatter3( XYZ1020(:,1), XYZ1020(:,2), XYZ1020(:,3), 'r','filled');
    % plot the matching MEG locations to these 10-20 locations
    H2 = scatter3( sens_unit(CHAN,1), sens_unit(CHAN,2), sens_unit(CHAN,3), 'b','filled');
    legend([H1(1),H2(1)],'10-20','MEG');
    view(2); rotate3d
end

if mod,
    % Edit ctf to reflect only these channels
    ctf.setup.number_channels = length(CHAN1020);
    ctf.data = ctf.data(:,CHAN,:);
    ctf.sensor.location = ctf.sensor.location(:,CHAN);
    ctf.sensor.orientation = ctf.sensor.orientation(:,CHAN);
    label = cell(1,length(CHAN1020));
    for i = 1:length(CHAN1020)
        label{i} = [ctf.sensor.label{CHAN(i)} '_' CHAN1020{i}];
    end
    ctf.sensor.label = label;
end

return
