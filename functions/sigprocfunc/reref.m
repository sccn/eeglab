% reref() - convert common reference EEG data to some other common reference
%           or to average reference
% Usage:
%   >> Dataout = reref(data);  % convert all channels to average reference
%   >> [Dataout Chanlocs Wout Sout meandata] = reref(data, refchan, 'key', 'val');
%                              % convert data to new reference with options
% Inputs:
%   data - 2-D data matrix (chans,frames*epochs) 
%   refchan  - reference channel number(s). There are two possibilities:
%          1) [] - if data state (see 'refstate' below) is 'common', compute 
%             average reference, otherwise convert from average reference.
%          2) [1 <= num(s) <= size(data,1)]: re-reference to channel num(s). 
%             If the 'withref' method (below) is used, compute the previous 
%             common reference channel as well.
% 
% Optional inputs:
%   'elocs'      - Current data electrode location structure (e.g., EEG.chanlocs).
%   'icaweights' - ICA weight matrix. Note: If this is ICA weights*sphere, 
%                  then the 'icasphere' input below should be [] or identity.
%   'icasphere'  - ICA sphere matrix (if any)
%   'icachansind' - Indices of the channels used in ICA decomposition
%
%   'method'     - ['standard'|'withref'] Do not ('standard') or do ('withref') 
%                  include reference channel data in output {def: 'standard'}. 
%                  Note: Option 'withref' not possible when multiple ref channel 
%                  indices are given as argument to 'refstate' (below).
%   'refstate  ' - ['common'|'averef'|[indices]] Current reference condition,
%                  ('averef') = average reference; ('common' or 0) = common
%                  reference. [indices] designate the current reference channel 
%                  or channels if present in the data {default: 'common'}
%   'refloc'     - Reference channel location {'label' theta radius} cell array
%                  containing the name and polar coordinates of the channel.
%                  If only 3-D coordinates are available, include the reference 
%                  channel location as the last channel in the 'elocs' structure.
%   'exclude'    - [integer array] channel indices to exclude from re-referencing
%                  (e.g., event marker channels, etc.)
%   'keepref'    - ['on'|'off'] keep reference channel in output (only usable 
%                  when there are several references).
%
% Outputs:
%   Dataout     - Input data converted to the new reference
%   Chanlocs    - Updated channel locations structure
%   Wout        - ICA weight matrix (former icaweights*icasphere)
%                 converted to new data reference
%   Sout        - ICA sphere matrix converted to an identity matrix
%   ICAinds     - New indices of channels used in ICA decomposition
%   meandata    - (1,dataframes) means removed from each data point
%
% Notes: 1) The average reference calculation implements two methods 
%           (see www.egi.com/Technotes/AverageReference.pdf)
%            V'i = (Vi-Vref) - sum(Vi-Vref)/number_of_electrodes
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
%
% Authors: Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, La Jolla, 1999-2002 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.32  2006/10/02 11:36:56  arno
% fix excluding channels
%
% Revision 1.30  2006/05/04 10:33:54  arno
% fixing last changes when not all channels are used for ICA
%
% Revision 1.28  2004/11/05 18:29:16  arno
% fixing channel label problem
%
% Revision 1.27  2003/12/11 17:55:22  arno
% remove debug msg
%
% Revision 1.26  2003/10/14 17:13:36  arno
% *** empty log message ***
%
% Revision 1.25  2003/10/14 17:12:07  arno
% *** empty log message ***
%
% Revision 1.24  2003/10/14 17:11:25  arno
% *** empty log message ***
%
% Revision 1.23  2003/07/29 18:39:59  arno
% debuging empty channel location structure
%
% Revision 1.22  2003/07/29 16:58:18  arno
% debuging exclude
%
% Revision 1.21  2003/07/28 16:44:26  arno
% allowing to include current ref channel
%
% Revision 1.20  2003/07/27 00:59:49  arno
% debuging re-referencing
%
% Revision 1.19  2003/07/26 01:13:56  arno
% debuging reref to several channels
%
% Revision 1.18  2003/07/26 00:01:32  arno
% brand new function
%
% Revision 1.17  2003/07/02 01:07:41  arno
% debug input check
%
% Revision 1.16  2003/06/11 00:10:27  arno
% debug multiple references
%
% Revision 1.15  2002/11/15 03:00:17  arno
% same
%
% Revision 1.14  2002/11/15 02:59:34  arno
% header for web
%
% Revision 1.13  2002/11/15 01:44:44  scott
% can not -> cannot
%
% Revision 1.12  2002/11/15 01:41:04  arno
% header for web
%
% Revision 1.11  2002/11/14 18:14:02  arno
% updating elocs parameter filter
%
% Revision 1.10  2002/11/13 23:04:41  arno
% removing extra electrode in average ref
%
% Revision 1.9  2002/11/13 22:28:57  arno
% removing the average reference with original common reference
%
% Revision 1.8  2002/11/13 20:29:41  arno
% debugging
%
% Revision 1.7  2002/11/13 15:12:21  scott
% help msg
%
% Revision 1.6  2002/11/13 15:08:40  scott
% help msg
%
% Revision 1.5  2002/11/12 23:32:59  arno
% debugging old ref channel potential
%
% Revision 1.4  2002/11/12 23:22:47  arno
% header typo
%
% Revision 1.3  2002/11/12 19:08:02  arno
% debugging
%
% Revision 1.2  2002/11/12 18:43:31  arno
% debug
%
% Revision 1.1  2002/11/12 17:58:08  arno
% Initial revision
%
% Revision 1.7  2002/09/05 00:30:23  scott
% added meandata output -sm
%
% Revision 1.6  2002/08/21 02:08:19  arno
% nothing
%
% Revision 1.5  2002/08/21 02:03:32  arno
% debugging ica reref
%
% Revision 1.4  2002/08/21 00:21:51  arno
% debugging
%
% Revision 1.3  2002/04/11 18:37:33  scott
% revised help msg
%
% Revision 1.2  2002/04/11 18:02:03  arno
% computing average reference of components
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 12/16/99 Corrected denomiator on the suggestion of Ian Nimmo-Smith, Cambridge UK
% 01-25-02 reformated help & license -ad 

function [data, Elocs, W, S, icachansind, meandata] = reref(data, ref, varargin)

if nargin<1
  help reref
  return
end
if nargin < 2
    ref = [];
end;

% check inputs
% ------------
g = finputcheck(varargin, { 'icaweight'   'real'    []          [];
                            'icaweights'  'real'    []          [];
                            'icasphere'   'real'    []          [];
                            'icachansind' 'integer'    []       [];
                            'method'     'string'  { 'standard' 'withref' }  'standard';
                            'refstate'   { 'string' 'integer' } { { 'common' 'averef' } [1 size(data,1)] }     'common';
                            'exclude'    'integer' [1 size(data,1)]          [];
                            'refloc'     'cell'    []                        {};
                            'keepref'    'string'  {'on' 'off' }             'off';
                            'elocs'      {'integer' 'struct'}  []            [] });
if isstr(g), error(g); end;
if ~isempty(g.icaweight)
    g.icaweights = g.icaweight;
end;
if ~isempty(g.icaweights) & ~isempty(g.icasphere)
    g.icaweights = g.icaweights*g.icasphere;
    % g.icasphere = eye(size(g.icasphere,1)); % make icasphere an identity matrix - see below
end;
if ~isempty(g.icaweights)
    if isempty(g.icachansind), 
        g.icachansind = [1:size(g.icaweights,2)]; 
        disp('Warning: reref() output has changed slightly since EEGLAB 5.02');
        disp('         the 4th output argument is the indices of channels used for ICA instead');
        disp('         of the mean reference value (which is now output argument 5)');
    end;
end;

chans = size(data,1);
if chans < 2 
  help reref
  return
end
if ref > chans+1
    error('reference channel index out of range');
end;

[dim1 dim2 dim3] = size(data);
data = reshape(data, dim1, dim2*dim3);

if ~isempty(g.exclude)
    rerefchans = setdiff([1:dim1], g.exclude);
    nbchans    = chans    - length(g.exclude);
else
    rerefchans = [1:dim1];
    nbchans    = chans;  
end;
rerefchansout = rerefchans;   

% return mean data
% ----------------
if nargout > 4
    meandata = sum(data(rerefchans,2))/nbchans;
end;

% compute average reference matrix
% --------------------------------
if ~isstr(g.refstate) | ~strcmp(g.refstate, 'averef')
    
    if strcmpi(g.method, 'withref')
        if ~isstr(g.refstate) & length(g.refstate) > 1
            error('Cannot include old reference as a channel when multiple reference channels are use');
        end;
        disp('Note: old reference electrode included in average reference computation');
        avematrix          = eye(nbchans)-ones(nbchans)*1/(nbchans+1); % reference is a relevant channel i.e. Cz
        avematrix(end+1,:) = -1/(nbchans+1);                           % potential for the new (previously ref) channel added to end of data
        if ~( length(g.elocs) > chans )
            if ~isempty(g.refloc) & ~isempty(g.refloc{1})
                if ~isempty(g.elocs)
                    g.elocs(end+1).labels = g.refloc{1};
                    g.elocs(end  ).theta  = g.refloc{2};
                    g.elocs(end  ).radius = g.refloc{3};
                else
                    disp('Reference-channel location ignored (no channel location structure provided)');
                end;
            else
                if ~isempty(g.elocs)
                    error('Need a old reference location to include it as a data channel');
                end;
            end;
        end;
        rerefchansout(end+1) = chans+1;
    else
        disp('Note: Unless present as a data channel, data for the old reference electrode are not included in the average reference computation');
        if length(g.refstate) > 1
            avematrix = eye(nbchans)-ones(nbchans)*1/nbchans;
        else
            avematrix = eye(nbchans)-ones(nbchans)*1/(nbchans-length(g.refstate)+1); % case of multiple references
        end;            
        if length(g.elocs) > chans
            g.elocs(end) = []; % discard info
        end;
    end;
else
    avematrix = eye(nbchans);
end;

% new reference electrode
% -----------------------
if strcmpi(g.method, 'withref')
    chans = chans+1;
end;

% generate rereferencing matrix
% -------------------------------
refmatrix = eye(chans); % begin with identiy matrix
if ~isempty(ref)    
    fprintf('Re-referencing data\n');
    for index = 1:length(ref)
        refmatrix(:,ref(index)) = refmatrix(:,ref(index))-1/length(ref);
    end;
end;

% deal with channel locations
% ------------------------------
if ~isempty(ref)    
    if length(ref) > 1 
        if strcmpi(g.keepref, 'off')
            % remove reference from channel location structure
            % -------------------------------------------------
            fprintf('Warning: reference channels have been removed\n');
            if ~isempty(g.elocs)
                g.elocs(ref) = [];
            end;
        end;
    else
        % copy channel location for single ref
        % ------------------------------------
        if ~isempty(g.elocs) & length(g.elocs) > chans & length(ref) == 1
            tmpelocs       = g.elocs(ref);
            g.elocs(ref)   = [];
            g.elocs(end+1) = tmpelocs;
        end;
    end;
end;

% exclude electrodes
% ------------------
rerefchansout = setdiff(rerefchansout, ref);
if strcmpi(g.keepref, 'on') & length(ref) > 1
    refmatrix(g.exclude,:) = []; % suppress excluded non-EEG channels
    refmatrix(:,g.exclude) = []; % suppress excluded non-EEG channels
else
    refmatrix([ref g.exclude],:) = []; % suppress references and non-EEG channels
    refmatrix(:,g.exclude      ) = []; % suppress non-EEG channels
end;

% rereferencing, there are faster methods
% but this one is the simpliest (for ICA re-referencing too)
% ------------------------------------------------------
data(rerefchansout,:) = (refmatrix*avematrix)*data(rerefchans,:); % perform matrix multiply
Elocs = g.elocs;
if strcmpi(g.keepref, 'off')
    rmchans = setdiff(1:size(data,1), ref);
    data  = data(rmchans,:);
end;

data = reshape(data, size(data,1), dim2, dim3);

% treat optional ica parameters
% -----------------------------
if ~isempty(g.icaweights) 

    % compute new channel indices
    % ---------------------------
    icachansind = g.icachansind;
    rmchans     = sort(rmchans);
    for i=length(icachansind):-1:1
        indchan = find( icachansind(i) == rmchans );
        if ~isempty( indchan )
            icachansind(i) = indchan;
        else
            icachansind(i) = [];
        end;
    end;

    % compute new weight matrix
    % -------------------------
    winv = pinv(g.icaweights);
    try, 
        W = pinv(refmatrix(icachansind,:)*avematrix(:,g.icachansind)*winv);
	catch,
        error([ 'Weight matrix size is different from the data size, re-referencing impossible' 10 ...
                      '(referencing must use the same number of channels (less excluded ones)' 10 ...
                      'as in the ICA decomposition. Best solution is to empty ICA weights, rereference, then run ICA decomposition again)']);
    end;
    S = eye(length(icachansind)); % output an identity sphere matrix
else
    W = []; S = [];
end;
