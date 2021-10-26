% loc_subsets() - Separate channels into maximally evenly-spaced subsets. 
%                 This is achieved by exchanging channels between subsets so as to
%                 increase the sum of average of distances within each channel subset.
% Usage:
%       >> subset = loc_subsets(chanlocs, nchans); % select an evenly spaced nchans
%       >> [subsets subidx pos] = loc_subsets(chanlocs, nchans, plotobj, plotchans, keepchans);
%
% Inputs:
%
%   chanlocs    - EEGLAB dataset channel locations structure (e.g., EEG.chanlocs)
%
%   nchans      - (1,N) matrix containing number of channels that should be in each 
%                 of N channel subsets respectively. if total number of channels in 
%                 all subsets is less than the number of EEG channels in the chanlocs 
%                 structure, N+1 subsets are created. The last subset includes the 
%                 remaining channels.
%
% Optional inputs:
%
%   plotobj     - (true|false|1|0) plot the time course of the objective function 
%                 {default: false|0) 
%   plotchans   - (true|false|1|0) make 3-D plots showing the channel locations of 
%                 the subsets {default: false|0) 
%   keepchans   - (cell array) channels that has to be kept in each set and
%                 not undergo optimization. You can use this option to make
%                 sure certain channels will be assigned to each set. For
%                 example, to keep channels 1:10 to the first subset and
%                 channels 20:30 to the second, use keepchans = {[1:10], [20:30]}.
%                 To only keeps channels 1:5 in the first set, use
%                 keepchans = {1:5}.  {default: {}}
%
% Outputs:
%
%   subsets     - {1,N} or {1,N+1} cell array containing channel indices for each 
%                 requested channel subset. 
%   subidx      - (1, EEG.nbchans) vector giving the index of the subset associated 
%                 with each channel.
%   pos         - (3,N) matrix, columns containing the cartesian (x,y,z) channel 
%                 positions plotted.
% Example:
%
%   % Create three sub-montages of a 256-channel montage containing 32, 60, and 100 
%   % channels respectively.The 64 remaining channels will be put into a fourth subset.  
%   % Also visualize time course of optimization and the channel subset head locations.
%
%   >> subset = loc_subsets(EEG.chanlocs, [32 60 100], true, true);
%
% Author: Nima Bigdely Shamlo, SCCN/INC/UCSD, 2007

% Copyright (C) 2007 Nima Bigdely Shamlo, SCCN/INC/UCSD, nima@sccn.ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% if plotOptimization|plotSubsets in line 130 removed by nima 3/6/2007
% line 133, removed num2str  removed by nima 3/6/2007

function [subset idx pos] = loc_subsets(chanlocs, numberOfChannelsInSubset, plotOptimization, plotSubsets, mandatoryChannelsForSet);

if sum(numberOfChannelsInSubset)> length(chanlocs)
    error('Total channels in requested subsets larger than number of EEG channels.');
end
if min(numberOfChannelsInSubset) < 2
    error('Number of channels in the requested subsets must be >= 2.');
end

rand('state',0);
if nargin < 5
    mandatoryChannelsForSet = {};
end

if nargin < 4
    plotSubsets = false;
end

if nargin < 3
    plotOptimization = false;
end

pos=[cell2mat({chanlocs.X}); cell2mat({chanlocs.Y}); cell2mat({chanlocs.Z});];
dist = squareform(pdist(pos'));

nChans = length(chanlocs);
idx = ones(nChans,1);
setId = cell(1, length(numberOfChannelsInSubset)); % cell array containing channels in each set

remainingChannels = 1:nChans; % channels that are to be assigned to subsets

% channels that have to stay in their original subset (as in
% mandatoryChannelsForSet) and should not be re-assigned
allMandatoryChannels = cell2mat(mandatoryChannelsForSet);

% assign requested mandatory channels to subset
for i=1:length(mandatoryChannelsForSet)
    setId{i} = mandatoryChannelsForSet{i};
    remainingChannels(mandatoryChannelsForSet{i}) = NaN; % flag with Nan so they can be deleted later, this is to keep indexing simple
end
remainingChannels(isnan(remainingChannels)) = [];

r = remainingChannels(randperm(length(remainingChannels)));

% randomly assign remaining channels to subsets
for i=1:length(numberOfChannelsInSubset)
    numberOfChannelsTobeAddedToSubset = numberOfChannelsInSubset(i) - length(setId{i})
    setId{i} = [setId{i} r(1:numberOfChannelsTobeAddedToSubset)];
    r(1:numberOfChannelsTobeAddedToSubset) = [];
end

if length(r) > 0
    setId{length(numberOfChannelsInSubset) + 1} = r; % last set gets remaining channels
end

if plotOptimization|plotSubsets
fprintf(['Creating total of ' num2str(length(setId)) ' channel subsets:\n']);
end

if plotOptimization
    figure;
  xp = floor(sqrt(length(setId)));
  yp = ceil(length(setId)/xp);
end
    counter = 1;
exchangeHappened = true;


while exchangeHappened
    exchangeHappened = false;
    for set1 = 1:length(setId)
        for set2 = (set1 + 1):length(setId)
            for i = 1:length(setId{set1})
                for j = 1:length(setId{set2})
                    chan(1) = setId{set1}(i);
                    chan(2) = setId{set2}(j);
                    if cost_of_exchanging_channels(chan,[set1 set2], setId, dist) < 0 && ~any(ismember(chan, allMandatoryChannels))
                        setId{set1}(find(setId{set1} == chan(1))) = chan(2);
                        setId{set2}(find(setId{set2} == chan(2))) = chan(1);
                        sumDistances(counter) = 0;
                        for s = 1:length(setId)
                            sumDistances(counter) = sumDistances(counter) ...
                                 + (sum(sum(dist(setId{s},setId{s}))) / length(setId{s}));
                        end
                        if plotOptimization
                            plot(1:counter,sumDistances,'-b');
                            xlabel('number of exchanges');
                            ylabel('sum mean distances within each channel subset');
                            drawnow;
                        else
                            if mod(counter, 20) ==0                             
                                fprintf('number of exchanges = %d\nsum of mean distances = %g\n',...
                                           counter, sumDistances(counter));
                            end
                        end
                        counter = counter + 1;
                        exchangeHappened = true;
                    end
                end
            end
        end
    end
end

for set = 1:length(setId)
    idx(setId{set}) = set;
%    legendTitle{set} = ['subset ' num2str(set)];
end

subset = setId;

if plotSubsets
    FIG_OFFSET = 40;
    fpos = get(gcf,'position');
    figure('position',[fpos(1)+FIG_OFFSET,fpos(2)-FIG_OFFSET,fpos(3),fpos(4)]);
    scatter3(pos(1,:), pos(2,:), pos(3,:),100,idx,'fill');
    axis equal;
    %legend(legendTitle); it does not work propr
    th=title('Channel Subsets');
    %set(th,'fontsize',14)
end

if length(r)>0 && (plotOptimization|plotSubsets)
    fprintf('The last subset returned contains the %d unused channels.\n',...
                   length(r));
end

function cost = cost_of_exchanging_channels(chan, betweenSets, setId, dist);

mirrorChan(1) = 2;
mirrorChan(2) = 1;

for i=1:2
    newSetId{betweenSets(i)} = setId{betweenSets(i)};
    newSetId{betweenSets(i)}(find(newSetId{betweenSets(i)} == chan(i))) = chan(mirrorChan(i));
end

cost = 0;
for i=betweenSets
    distSumBefore{i} = sum(sum(dist(setId{i},setId{i})));
    distSumAfter{i} = sum(sum(dist(newSetId{i},newSetId{i})));
    cost = cost + (distSumBefore{i} - distSumAfter{i}) / length(setId{i});
end


%cost = (distSumAfter{1} > distSumBefore{1}) && (distSumAfter{2} > distSumBefore{2});

