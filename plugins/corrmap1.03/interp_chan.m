% interp_chan() - interpolates missing channels in datasets stored in
%                 ALLEEG based on channel information provided and returns
%                 cell array containg EEG.icawinv with the same
%                 dimensions for these datasets
%
% Usage:
%             >> comp_interp=interp_chan(ALLEEG,chanlocs);
% Inputs:
%   ALLEEG    - structure containing all datasets
%   chanlocs  - path for a channel locations file
%
% Outputs:
%   comp_interp   - a cell array containg EEG.icawinv interpolated for all
%   datasets
%
%  See also:  corrmap(), pop_corrmap() and eeg_interp()
%
% Authors: Filipa Campos Viola, 25/01/2008, MRC-IHR, Southampton, UK
% (f.viola@soton.ac.uk)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) F. Campos Viola, MRC-IHR
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

% revised by F Campos-Viola - corrmap1.01 (30/01/2009)

%Updated - 29/01/2008 (updated - method was changed from 'default' to 'v4=invdist'

function [comp_interp,CHANS]=interp_chan(ALLEEG,chanlocs)

CHANS=[];
%chanlocs='Q:\\StefansData\\asym01\\config\\68_channel_besa_sphere.elp';

%reading good channels from channel location file
CHANS=pop_chanedit(CHANS,'load',{chanlocs, 'filetype', 'autodetect'});

aux=1:length(ALLEEG);

comp={ALLEEG(aux).icawinv};

comp_interp=cell(1,length(comp)); %cell to store new comp

x_comp=zeros(1,length(CHANS));
y_comp=zeros(1,length(CHANS));

for i=1:length(comp)

    temp=zeros(length(CHANS),1);

    [xgood,ygood] = pol2cart([CHANS.theta],[CHANS.radius]); %positions for all channels

    [x_orig,y_orig] = pol2cart([ALLEEG(i).chanlocs.theta],[ALLEEG(i).chanlocs.radius]);
    %positions for channels recorded in this dataset


    for j=1:size(ALLEEG(i).icawinv,2)

        if ~isempty(ALLEEG(i).badchannels)
            %updating variable "temp" where interpolated data will be stored
            h=1;
            g=1;
            while h<=length(CHANS)
                aux=find(ALLEEG(i).badchannels==h);
                if isempty(aux)
                    temp(h,1)=ALLEEG(i).icawinv(g,j);
                    h=h+1;
                    g=g+1;
                else
                    temp(h,1)=0;
                    h=h+1;
                    g=g;
                end

            end

            aaa = griddata(y_orig,x_orig,ALLEEG(i).icawinv(:,j),...
                ygood(ALLEEG(i).badchannels),xgood(ALLEEG(i).badchannels),'v4'); % interpolate data

            clear x_comp
            clear y_comp

            temp(ALLEEG(i).badchannels)=aaa';

            comp_interp{i}(:,j)=temp;
            clear temp

        else
            comp_interp{i}(:,j)=ALLEEG(i).icawinv(:,j); % bad channes don't exist, so no need to interpolate
        end

    end

end

