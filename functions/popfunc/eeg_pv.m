% eeg_pv()   - Compute EEG.data 'percent variance ' (pv) of the whole EEG data versus the projections
%              of specified components. 
%              Can omit specified components and channels from the computation. Can draw a plot 
%              of the scalp distribution of pv, or progressively compute the pv for comps
%              1:k, where k = 1 -> the total number of components.  Note: pv's of spatially
%              non-orthogonal independent components may not add to 100%, and individual component 
%              pv could be < 0%.
% Usage:
%              >> [pv] = eeg_pv(EEG,comps);
%              >> [pv,pvs,vars] = eeg_pv(EEG,comps,artcomps,omitchans,fraction,'plot');
% Inputs:
%    EEG       - EEGLAB dataset. Must have icaweights, icasphere, icawinv, icaact.
%    comps     - vector of component indices to sum {default|[] -> progressive mode}
%                In progressive mode, comps is first [1], then [1 2], etc. up to
%                [1:size(EEG.icaweights,2)] (all components); here, the plot shows pv.
%    artcomps  - vector of artifact component indices to remove from data before
%                computing pv {default|[]: none}
%    omitchans - channels to omit from the computation (e.g. off-head, etc) {default|[]: none}
%    fraction  - (0<real<=1) fraction of the data to randomly select {default|[]: 1=all}
%    'plot'    - Plot scalp map of channel pvs. {default: Plot only if no output arguments}
%
% Outputs:
%    pv      - (real) percent total variance accounted for by the summed back-projection of
%                the requested components. If comps is [], a vector of pvs for the sum of 
%                components 1:k (k=1:ncomps).
%    pvs     - (real vector) percent variance accounted for by the summed back-projection of
%                the requested components to each data channel. If comps is [], a matrix of
%                pvs (as for pv above).
%    vars      - variances of the requested channels
%
% Fields:  
%    Assumes existence of the following EEG fields: EEG.data, EEG.pnts, EEG.nbchan, EEG.trials,
%          EEG.icaact, EEG.icaweights, EEG.icasphere, EEG.icawinv, and for plot, EEG.chanlocs
%
% See also:  eeg_pvaf()
%
% Author: from eeg_pvaf(), Scott Makeig, SCCN/INC/UCSD, 02/04/05

function [pv,pvs,pvall] = eeg_pv(EEG,comps,artcomps,omitchans,fraction,plotflag)

if nargin < 1 | nargin > 6
   help eeg_pv
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
  fprintf('eeg_pv(): considering all the data.\n');
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
  help eeg_pv
  error('at least one channel to omit > number of channels in data');
 end
 if min(omitchans)<1
  help eeg_pv
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
   help eeg_pv
   return
end
if isempty(EEG.icasphere)
   help eeg_pv
   return
end

if isempty(EEG.icawinv)
     EEG.icawinv = pinv(EEG.icaweights*EEG.icasphere);
end
if isempty(EEG.icaact)
    help eeg_pv
    fprintf('EEG.icaact not present.\n');
    % EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data; % remake it like this
end
if max(comps) > size(EEG.icawinv,1)
    help eeg_pv
   fprintf('Only %d components in this dataset. Cannot project component %d.\n',numcomps,max(comps));
   error('bad comps input');
end
if ~isempty(artcomps) & max(artcomps) > numcomps
    help eeg_pv
   fprintf('Only %d components in this dataset. Cannot project artcomp %d.\n',numcomps,max(artcomps));
   error('bad artcomps input')
end

npts = EEG.trials*EEG.pnts;
allcomps = 1:numcomps;
if progressive
   fprintf('Considering components up to: ');
   cum_pv  = zeros(1,numcomps);
   cum_pvs = zeros(numcomps,numchans);
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
          fprintf('eeg_pv(): not back-projecting %d comps already in the artcomps.\n',length(c));
        else
          fprintf('eeg_pv(): not back-projecting comp %d already in the artcomps.\n',comps(c));
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

chans
comps
size(EEG.icawinv(chans,comps))
size(EEG.icaact(comps,:)')
pvcomp =  var((EEG.icawinv(chans,comps)*EEG.icaact(comps,:))', varwts);

%
%%%%%%%%%%%%%%%%%%%%%%%% compute percent variance %%%%%%%%%%%%%%%
%
pvs =  pvcomp ./ pvall;
pvs = 100*pvs;
pv = sum(pvcomp) ./ sum(pvall);
pv = 100*pv;

if ~progressive
  break
else
  cum_pv(comp) = pv;
  cum_pvs(comp,:) = pvs;
end

end %%%%%%%%%%%%%% end progressive forloop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if progressive % output accumulated results
  fprintf('\n');
  pv = cum_pv;
  pvs = cum_pvs;
  if plotit
     plot(1:numcomps,pv);
     xl = xlabel('Components Included (1:n)');
     yl = ylabel('Percent Variance Accounted For (pv)');
     set(xl,'fontsize',15);
     set(yl,'fontsize',15);
     set(gca,'fontsize',14);
  end
elseif plotit 
%
%%%%%%%%%%%%%%%%%%%%%%%% plot the scalp distribtion of pv %%%%%%%%%%%%%
%
 if isfield(EEG,'chanlocs')
   chanlocs = EEG.chanlocs;
   if ~isempty(omitchans)
     chanlocs(omitchans) = [];
   end
   topoplot(pvs',chanlocs);  % plot pv here

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
   if max(pvs)>100, 
      maxc=max(pvs)
   else 
      maxc=100; 
   end;

   pvstr=sprintf('Total pv: %3.1f%%',pv);
   tx=text(-0.9,-0.6,pvstr);

   caxis([-100 100]);
   cb=cbar('vert',33:64,[0 100]); % color bar showing >0 (green->red) only
 else
    fprintf('EEG.chanlocs not found - not plotting scalp pv\n');
 end
end % plotit


