% reref() - convert common reference EEG data to some other common reference
%           or to average reference
% Usage:
%   >> Dataout = reref(data);  % convert all channels to average reference
%   >> [Dataout Chanlocs] = reref(data, refchan, 'key', 'val');
%                              % convert data to new reference with options
% Inputs:
%   data - 2-D or 3-D data matrix (chans,frames*epochs) 
%   refchan  - reference channel number(s). There are three possibilities:
%          1) [] -  compute average reference
%          2) [X]: re-reference to channel X 
%          2) [X Y Z ...]: re-reference to the average of channel X Y Z ... 
% 
% Optional inputs:
%   'exclude'    - [integer array] channel indices to exclude from re-referencing
%                  (e.g., event marker channels, etc.)
%   'keepref'    - ['on'|'off'] keep reference channel in output (only usable 
%                  when there are several references).
%   'elocs'      - Current data electrode location structure (e.g., EEG.chanlocs).
%   'refloc'     - Reference channel location single element structure or cell array
%                  {'label' theta radius} containing the name and polar coordinates 
%                  of the current channel. Including this entry means that
%                  this channel will be included (reconstructed) in the
%                  output.
%
% Outputs:
%   Dataout     - Input data converted to the new reference
%   Chanlocs    - Updated channel locations structure
%
% Notes: 1) The average reference calculation implements two methods 
%           (see www.egi.com/Technotes/AverageReference.pdf)
%            V'i = (Vi-Vref) - sum(Vi-Vref)/number_of_electrodes
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2009-
%          previous version: Arnaud Delorme & Scott Makeig, 1999-2002

% Deprecated inputs:
% These inputs are still accepted but not processed. The function returns
% accurate results irrespective of the entry of these options.
%   'refstate  ' - ['common'|'averef'|[indices]] Current reference condition,
%                  ('averef') = average reference; ('common' or 0) = common
%                  reference. [indices] designate the current reference channel 
%                  or channels if present in the data {default: 'common'}
%   'method'     - ['standard'|'withref'] Do not ('standard') or do ('withref') 
%                  include reference channel data in output {def: 'standard'}. 
%                  Note: Option 'withref' not possible when multiple ref channel 
%                  indices are given as argument to 'refstate' (below).
%
% ICA inputs:
% These inputs are still accepted but not the ICA conversion is now
% performed from within pop_reref()
%   'icaweights' - ICA weight matrix. Note: If this is ICA weights*sphere, 
%                  then the 'icasphere' input below should be [] or identity.
%   'icasphere'  - ICA sphere matrix (if any)
%   'icachansind' - Indices of the channels used in ICA decomposition
%
% Outputs:
%   Wout        - ICA weight matrix (former icaweights*icasphere)
%                 converted to new data reference
%   Sout        - ICA sphere matrix converted to an identity matrix
%   ICAinds     - New indices of channels used in ICA decomposition
%   meandata    - (1,dataframes) means removed from each data point
%
%        2) In conversion of the weight matrix to a new reference
%           where WS = Wts*Sph  and  ica_act = WS*data, then
%              data = inv(WS)*ica_act;
%           If R*data are the re-referenced data,
%              R*data= R*inv(WS)*ica_act;   
%           And Wout = inv(R*inv(WS));
%           Now, Sout = eye(length(ICAinds));
%           The re-referenced ICA component maps are now the 
%           columns of inv(Wout), and the icasphere matrix, Sout, 
%           is an identity matrix. Note: inv() -> pinv() when 
%           PCA dimension reduction is used during ICA decomposition.

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 12/16/99 Corrected denomiator on the suggestion of Ian Nimmo-Smith, Cambridge UK
% 01-25-02 reformated help & license -ad 

function [data, Elocs, morechans, W, S, icachansind, meandata] = reref(data, ref, varargin)

if nargin<1
  help reref
  return
end
if nargin < 2
    ref = [];
end

% check inputs
% ------------
g = finputcheck(varargin, { 'icaweight'   'real'    []          [];
                            'icaweights'  'real'    []          [];
                            'icasphere'   'real'    []          [];
                            'icachansind' 'integer'    []       [];
                            'interpchan'  {''}      []          [];
                            'method'     'string'  { 'standard','withref' }  'standard';
                            'refstate'   { 'string','integer' } { { 'common','averef' } [1 size(data,1)] }     'common'; % ot used but kept for backward compatib.
                            'exclude'    'integer' [1 size(data,1)]          [];
                            'refloc'     { 'cell','struct' }  { [] [] }   {};
                            'keepref'    'string'  {'on','off' }             'off';
                            'elocs'      {'integer','struct'}  []            [] });
if ischar(g), error(g); end
if ~isempty(g.icaweight)
    g.icaweights = g.icaweight;
end
if ~isempty(g.icaweights)
    if isempty(g.icachansind), 
        g.icachansind = [1:size(g.icaweights,2)]; 
        disp('Warning: reref() output has changed slightly since EEGLAB 5.02');
        disp('         the 4th output argument is the indices of channels used for ICA instead');
        disp('         of the mean reference value (which is now output argument 5)');
    end
end
if ischar(ref), ref = { ref }; end
if iscell(ref), ref = eeg_chaninds(g.elocs, ref); end
if ~isempty(ref)
    if ref > size(data,1)
        error('reference channel index out of range');
    end
end

[dim1 dim2 dim3] = size(data);
data = reshape(data, dim1, dim2*dim3);

% single reference not present in the data
% add it as blank data channel at the end
% ----------------------------------------
if ~isempty(g.refloc) == 1
    if ~isempty(g.elocs)
        if iscell(g.refloc)
            data(end+1,:) = 0;
            g.elocs(end+1).labels = g.refloc{1};
            g.elocs(end  ).theta  = g.refloc{2};
            g.elocs(end  ).radius = g.refloc{3};
        else
            data(end+length(g.refloc),:) = 0;
            for iLocs = 1:length(g.refloc)
                g.elocs(end+1).labels = g.refloc(iLocs).labels;
                fieldloc = fieldnames(g.elocs);
                for ind = 1:length(fieldloc)
                    g.elocs(end) = setfield(g.elocs(end), fieldloc{ind}, getfield(g.refloc(iLocs), fieldloc{ind}));
                end
            end
        end
    end
    [dim1 dim2 dim3] = size(data);
end

% exclude some channels
% ---------------------
chansin   = setdiff_bc([1:dim1], g.exclude);
nchansin  = length(chansin);

% return mean data
% ----------------
if nargout > 4
    meandata = sum(data(chansin,2))/nchansin;
end

% generate rereferencing matrix
% -----------------------------
if 0 % alternate code - should work exactly the same
    if isempty(ref)
        ref=chansin; % average reference
    end % if 
    chansout=chansin; 
    data(chansout,:)=data(chansout,:)-ones(nchansin,1)*mean(data(ref,:),1); 
else
    if ~isempty(ref) % not average reference   
        refmatrix = eye(nchansin); % begin with identity matrix
        tmpref = ref;
        for index = length(g.exclude):-1:1
            tmpref(find(g.exclude(index) < tmpref)) = tmpref(find(g.exclude(index) < tmpref))-1;
        end
        for index = 1:length(tmpref)
            refmatrix(:,tmpref(index)) = refmatrix(:,tmpref(index))-1/length(tmpref);
        end
    else % compute average reference
        refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin;
    end
    chansout = chansin;
    data(chansout,:) = refmatrix*data(chansin,:);
end

% change reference in elocs structure
% -----------------------------------
if ~isempty(g.elocs)
    if isempty(ref)
        for ind = chansin
            g.elocs(ind).ref = 'average';
        end
    else
        reftxt = { g.elocs(ref).labels };
        if length(reftxt) == 1, reftxt = reftxt{1}; 
        else
            reftxt = cellfun(@(x)([x ' ']), reftxt, 'uniformoutput', false);
            reftxt = [ reftxt{:} ];
        end
        for ind = chansin
            g.elocs(ind).ref = reftxt;
        end
    end
end

% remove reference
% ----------------
morechans = [];
if strcmpi(g.keepref, 'off')
    data(ref,:) = [];
    if ~isempty(g.elocs)
        morechans = g.elocs(ref);
        g.elocs(ref) = [];
    end
end

data = reshape(data, size(data,1), dim2, dim3);

% treat optional ica parameters
% -----------------------------
W = []; S = []; icachansind = [];
if ~isempty(g.icaweights) 
    disp('Warning: This function does not process ICA array anymore, use the pop_reref function instead');
end
Elocs = g.elocs;
