% eeg_pvaf() - Compute EEG.data 'percent variance accounted for' (pvaf) by specified components. 
%              Can omit specified components and channels from the computation. Can draw a plot 
%              of the scalp distribution of pvaf, or progressively compute the pvaf for comps
%              1:k, where k = 1 -> the total number of components.  Note: pvaf's of spatially
%              non-orthogonal independent components may not add to 100%, and individual component 
%              pvaf could be < 0%.
% Usage:
%              >> [pv] = eeg_pvaf(EEG,comps);s
%              >> [pvaf,pvafs,vars] = eeg_pvaf(EEG, comps,'key', val);
% Inputs:
%    EEG       - EEGLAB dataset. Must have icaweights, icasphere, icawinv, icaact.
%    comps     - vector of component indices to sum {default|[] -> progressive mode}
%                In progressive mode, comps is first [1], then [1 2], etc. up to
%                [1:size(EEG.icaweights,2)] (all components); here, the plot shows pvaf.
%
% Optional inputs:
%    'artcomps'  - [integer] vector of artifact component indices to remove from data before
%                  computing pvaf {default|[]: none}
%    'omitchans' - [integer] channels to omit from the computation (e.g. off-head, etc) 
%                  {default|[]: none}
%    'chans'     - [integer] only compute pvaf at selected channels. Overwrite omitchans above.
%    'fraction'  - [0<real<=1] fraction of the data to randomly select {default|[]: 1=all}
%    'plot'      - ['on'|'off'] Plot scalp map of channel pvafs. {default: Plot only if no 
%                  output arguments}
%
% Outputs:
%    pvaf      - (real) percent total variance accounted for by the summed back-projection of
%                the requested components. If comps is [], a vector of pvafs for the sum of 
%                components 1:k (k=1:ncomps).
%    pvafs     - (real vector) percent variance accounted for by the summed back-projection of
%                the requested components to each data channel. If comps is [], a matrix of
%                pvafs (as for pv above).
%    vars      - variances of the requested channels
%
% Fields:  
%    Assumes existence of the following EEG fiels: EEG.data, EEG.pnts, EEG.nbchan, EEG.trials,
%          EEG.icaact, EEG.icaweights, EEG.icasphere, EEG.icawinv, and for plot, EEG.chanlocs
%
% Author: Scott Makeig & Arnaud Delorme, SCCN, INC, UCSD, Fri Feb 13, 2004

% Copyright (C) 2004- Scott Makeig & Arnaud Delorme, SCCN, INC, UCSD
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

function [pvaf,pvafs,pvall] = eeg_pvaf(EEG,comps, varargin)

if nargin < 1
   help eeg_pvaf
   return
end

g = finputcheck(varargin, { 'artcomps'   'integer'    []         [];
                            'omitchans'  'integer'    []         [];
                            'chans'      'integer'    []         [];
                            'fraction'   'real'       []         1;
                            'plot'       'string'     { 'on';'off';'def' } 'def' }, 'eeg_pvaf');
if isstr(g), error(g); end;

numcomps = size(EEG.icaact,1);
if round(g.fraction*EEG.pnts*EEG.trials)<1
   error('g.fraction of data specified too small.')
   return
end
if strcmpi(g.plot, 'def')
    if nargout > 0, g.plot = 'on';
    else            g.plot = 'off';
    end;
end 

numchans = EEG.nbchan;
chans = 1:numchans;
if ~isempty(g.chans)
    g.omitchans = setdiff([1:EEG.nbchan], g.chans);
end;
if ~isempty(g.omitchans)
 if max(g.omitchans)>numchans
  help eeg_pvaf
  error('at least one channel to omit > number of channels in data');
 end
 if min(g.omitchans)<1
  help eeg_pvaf
  error('channel numbers to omit must be > 0');
 end
 chans(g.omitchans) = [];
end

progressive = 0; % by default, progressive mode is off
if nargin < 2 | isempty(comps)|comps==0
  comps = [];
  progressive = 1;  % turn progressive mode on
end

if isempty(EEG.icaweights)
   help eeg_pvaf
   return
end
if isempty(EEG.icasphere)
   help eeg_pvaf
   return
end

if isempty(EEG.icawinv)
     EEG.icawinv = pinv(EEG.icaweights*EEG.icasphere);
end
if isempty(EEG.icaact)
    help eeg_pvaf
    fprintf('EEG.icaact not present.\n');
    % EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data; % remake it like this
end
if max(comps) > size(EEG.icawinv,1)
    help eeg_pvaf
   fprintf('Only %d components in this dataset. Cannot project component %d.\n',numcomps,max(comps));
   error('bad comps input');
end
if ~isempty(g.artcomps) & max(g.artcomps) > numcomps
    help eeg_pvaf
   fprintf('Only %d components in this dataset. Cannot project artcomp %d.\n',numcomps,max(g.artcomps));
   error('bad artcomps input')
end

npts = EEG.trials*EEG.pnts;
allcomps = 1:numcomps;
if progressive
   fprintf('Considering components up to: ');
   cum_pvaf  = zeros(1,numcomps);
   cum_pvafs = zeros(numcomps,numchans);
end

for comp = 1:numcomps %%%%%%%%%%%%%%% progressive mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if progressive
   comps = allcomps(1:comp); % summing components 1 to current comp
   fprintf('%d ',comp)
end

if ~isempty(g.artcomps)
   [a b c] = intersect(g.artcomps,comps);
   if ~isempty(a)
      if ~progressive
        if length(a)>1
          fprintf('eeg_pvaf(): not back-projecting %d comps already in the artcomps.\n',length(c));
        else
          fprintf('eeg_pvaf(): not back-projecting comp %d already in the artcomps.\n',comps(c));
        end
      end
      comps(c) = [];
   end
end
if ~isempty(g.artcomps) & min([comps g.artcomps]) < 1
   error('comps and artcomps must contain component indices');
end

%
%%%%%%%%%%%%%%%%%%%%%%%% compute variance accounted for by specified components %%%%%%%%%%%%%
%
if ~progressive | comp == 1 % pare out g.omitchans and artcomps from EEG.data
  if ~isempty(g.artcomps)
    EEG.data = EEG.data(chans,:) - EEG.icawinv(chans,g.artcomps)*EEG.icaact(g.artcomps,:);
  else
    EEG.data = EEG.data(chans,:); 
  end

  nsel = round(g.fraction*npts);
  varpts = randperm(npts);
  varwts = ones(size(varpts));
  if nsel<npts
    varwts(varpts(nsel+1:npts)) = 0;
  end
  pvall = var(EEG.data(:,:)',varwts);
end

pvdiff =  var((EEG.data(:,:) - EEG.icawinv(chans,comps)*EEG.icaact(comps,:))', varwts);

%
%%%%%%%%%%%%%%%%%%%%%%%% compute percent variance accounted for %%%%%%%%%%%%%%%
%
pvafs =  pvdiff ./ pvall;
pvafs = 100-100*pvafs;
pvaf = sum(pvdiff) ./ sum(pvall);
pvaf = 100-100*pvaf;

if ~progressive
  break
else
  cum_pvaf(comp) = pvaf;
  cum_pvafs(comp,:) = pvafs;
end

end %%%%%%%%%%%%%% end progressive forloop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if progressive % output accumulated results
  fprintf('\n');
  pvaf = cum_pvaf;
  pvafs = cum_pvafs;
  if strcmpi(g.plot, 'on');
     plot(1:numcomps,pvaf);
     xl = xlabel('Components Included (1:n)');
     yl = ylabel('Percent Variance Accounted For (pvaf)');
     set(xl,'fontsize',15);
     set(yl,'fontsize',15);
     set(gca,'fontsize',14);
  end
elseif strcmpi(g.plot, 'on')
%
%%%%%%%%%%%%%%%%%%%%%%%% plot the scalp distribtion of pvaf %%%%%%%%%%%%%
%
 if isfield(EEG,'chanlocs')
   chanlocs = EEG.chanlocs;
   if ~isempty(g.omitchans)
     chanlocs(g.omitchans) = [];
   end
   if length(chanlocs) > 1
       topoplot(pvafs',chanlocs);  % plot pvaf here
   end;
   
   if length(comps)>5        % add text legend
     if length(g.artcomps)>3
        tlstr=sprintf('Pvaf by %d comps in data minus %d comps',length(comps),length(g.artcomps));
     elseif isempty(g.artcomps)
        tlstr=sprintf('Pvaf by %d comps in data',length(comps));
     elseif length(g.artcomps)==1 % < 4 g.artcomps, list them
        tlstr=sprintf('Pvaf by %d comps in data (less comp ',length(comps));
        tlstr = [tlstr sprintf('%d ',g.artcomps) ')'];
     else
        tlstr=sprintf('Pvaf by %d comps in data (less comps ',length(comps));
        tlstr = [tlstr sprintf('%d ',g.artcomps) ')'];
     end
   else %  < 6 comps, list them
     if length(comps)>1
        tlstr=sprintf('Pvaf by comps ');
     else
        tlstr=sprintf('Pvaf by comp ');
     end
     if length(g.artcomps)>3
        tlstr = ...
[tlstr sprintf('%d ',comps) sprintf('in data minus %d comps',length(comps),length(g.artcomps))];
     else
        if isempty(g.artcomps)
           tlstr = [tlstr sprintf('%d ',comps) 'in data'];
        elseif length(g.artcomps)==1
           tlstr = [tlstr sprintf('%d ',comps) 'in data (less comp '];
           tlstr = [tlstr int2str(g.artcomps) ')'];
        else
           tlstr = [tlstr sprintf('%d ',comps) 'in data (less comps '];
           tlstr = [tlstr sprintf('%d ',g.artcomps) ')'];
        end
     end
   end
   tl=title(tlstr);
   if max(pvafs)>100, 
      maxc=max(pvafs)
   else 
      maxc=100; 
   end;

   pvstr=sprintf('Total pvaf: %3.1f%%',pvaf);
   tx=text(-0.9,-0.6,pvstr);

   caxis([-100 100]);
   cb=cbar('vert',33:64,[0 100]); % color bar showing >0 (green->red) only
 else
    fprintf('EEG.chanlocs not found - not plotting scalp pvaf\n');
 end
end % end plot


