% eeg_pvaf() - Compute EEG.data 'percent variance accounted for' (pvaf) by specified components. 
%              Can omit specified components and channels from the computation. Can draw a plot 
%              of the scalp distribution of pvaf, or progressively compute the pvaf for comps
%              1:k, where k = 1 -> the total number of components.  Note: pvaf's of spatially
%              non-orthogonal independent components may not add to 100%, and individual component 
%              pvaf could be < 0%.
% Usage:
%              >> [pv] = eeg_pvaf(EEG,comps);
%              >> [pvaf,pvafs,vars] = eeg_pvaf(EEG,comps,artcomps,omitchans,fraction,'plot');
% Inputs:
%    EEG       - EEGLAB dataset. Must have icaweights, icasphere, icawinv, icaact.
%    comps     - vector of component indices to sum {default|[] -> progressive mode}
%                In progressive mode, comps is first [1], then [1 2], etc. up to
%                [1:size(EEG.icaweights,2)] (all components); here, the plot shows pvaf.
%    artcomps  - vector of artifact component indices to remove from data before
%                computing pvaf {default|[]: none}
%    omitchans - channels to omit from the computation (e.g. off-head, etc) {default|[]: none}
%    fraction  - (0<real<=1) fraction of the data to randomly select {default|[]: 1=all}
%    'plot'    - Plot scalp map of channel pvafs. {default: Plot only if no output arguments}
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

% Author: Scott Makeig, SCCN/INC/UCSD, Fri Feb 13, 2004

function [pvaf,pvafs,pvall] = eeg_pvaf(EEG,comps,artcomps,omitchans,fraction,plotflag)

if nargin < 1 | nargin > 6
   help eeg_pvaf
   return
end
numcomps = size(EEG.icaact,1);
plotit = 0;
if nargin>5 | nargout < 1
   plotit = 1;
end
if nargin<5 | isempty(fraction)
  fraction = 1;
end
if fraction>1
  fprintf('eeg_pvaf(): considering all the data.\n');
  fraction=1;
end
if round(fraction*EEG.pnts*EEG.trials)<1
   error('fraction of data specified too small.')
   return
end
if nargin<4 | isempty(omitchans)
  omitchans = [];
end
if nargin<3|isempty(artcomps)
  artcomps=[];
end

numchans = EEG.nbchan;
chans = 1:numchans;
if ~isempty(omitchans)
 if max(omitchans)>numchans
  help eeg_pvaf
  error('at least one channel to omit > number of channels in data');
 end
 if min(omitchans)<1
  help eeg_pvaf
  error('channel numbers to omit must be > 0');
 end
 chans(omitchans) = [];
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
if ~isempty(artcomps) & max(artcomps) > numcomps
    help eeg_pvaf
   fprintf('Only %d components in this dataset. Cannot project artcomp %d.\n',numcomps,max(artcomps));
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

if ~isempty(artcomps)
   [a b c] = intersect(artcomps,comps);
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
if ~isempty(artcomps) & min([comps artcomps]) < 1
   error('comps and artcomps must contain component indices');
end

%
%%%%%%%%%%%%%%%%%%%%%%%% compute variance accounted for by specified components %%%%%%%%%%%%%
%
if ~progressive | comp == 1 % pare out omitchans and artcomps from EEG.data
  if ~isempty(artcomps)
    EEG.data = EEG.data(chans,:) - EEG.icawinv(chans,artcomps)*EEG.icaact(artcomps,:);
  else
    EEG.data = EEG.data(chans,:); 
  end

  nsel = round(fraction*npts);
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
  if plotit
     plot(1:numcomps,pvaf);
     xl = xlabel('Components Included (1:n)');
     yl = ylabel('Percent Variance Accounted For (pvaf)');
     set(xl,'fontsize',15);
     set(yl,'fontsize',15);
     set(gca,'fontsize',14);
  end
elseif plotit 
%
%%%%%%%%%%%%%%%%%%%%%%%% plot the scalp distribtion of pvaf %%%%%%%%%%%%%
%
 if isfield(EEG,'chanlocs')
   chanlocs = EEG.chanlocs;
   if ~isempty(omitchans)
     chanlocs(omitchans) = [];
   end
   topoplot(pvafs',chanlocs);  % plot pvaf here

   if length(comps)>5        % add text legend
     if length(artcomps)>3
        tlstr=sprintf('Pvaf by %d comps in data minus %d comps',length(comps),length(artcomps));
     elseif isempty(artcomps)
        tlstr=sprintf('Pvaf by %d comps in data',length(comps));
     elseif length(artcomps)==1 % < 4 artcomps, list them
        tlstr=sprintf('Pvaf by %d comps in data (less comp ',length(comps));
        tlstr = [tlstr sprintf('%d ',artcomps) ')'];
     else
        tlstr=sprintf('Pvaf by %d comps in data (less comps ',length(comps));
        tlstr = [tlstr sprintf('%d ',artcomps) ')'];
     end
   else %  < 6 comps, list them
     if length(comps)>1
        tlstr=sprintf('Pvaf by comps ');
     else
        tlstr=sprintf('Pvaf by comp ');
     end
     if length(artcomps)>3
        tlstr = ...
[tlstr sprintf('%d ',comps) sprintf('in data minus %d comps',length(comps),length(artcomps))];
     else
        if isempty(artcomps)
           tlstr = [tlstr sprintf('%d ',comps) 'in data'];
        elseif length(artcomps)==1
           tlstr = [tlstr sprintf('%d ',comps) 'in data (less comp '];
           tlstr = [tlstr int2str(artcomps) ')'];
        else
           tlstr = [tlstr sprintf('%d ',comps) 'in data (less comps '];
           tlstr = [tlstr sprintf('%d ',artcomps) ')'];
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
end % plotit


