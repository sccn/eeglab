% reref() - convert common reference EEG data to some other common reference
%           or to average reference
%
% Usage:
%   >> dataout = reref(data);
%   >> [dataout Chanlocs Wout Sout meandata] = reref(data, refchan, 'key', 'val');
%
% Inputs:
%   data - 2-D data matrix (chans,frames*epochs) 
%   refchan  - reference channel number(s) -- two possibilities here
%          1) [] - if data state (see 'refstate' below) is 'common' compute 
%             average reference, otherwise undo average reference.
%          2) 1 <= num <= size(data,1): re-reference to channel num. If the 
%             'withref' method is set, the function computes the previous 
%             common reference channel.
% 
% Optional inputs:
%   'elocs'      - Electrode location structure (e.g., EEG.chanlocs).
%   'icaweight'  - ICA weight matrix (Note: this should be weights*sphere)
%   'method'     - ['standard'|'withref'] Include ('withref') reference channel 
%                  in output. Default is 'standard'. 
%   'refstate  ' - ['common'|'averef'|integer] Current average reference.
%                  Use this parameter to re-reference data to a given channel if data
%                  is already average referenced (by setting it to 'averef').
%                  Integer designates channel reference indices.
%                  Default is 'common'. 
%   'refloc'     - Reference channel location -- a cell array with name and 
%                  polar coordinates of the channel, { 'label' theta radius }. 
%                  For 3-D location, include the reference as the last channel 
%                  in the 'chanlocs' structure.
%   'exclude'    - [integer array] channel indices to exclude from rereferencing
%                  i.e. EMG
%   'keepref'    - ['on'|'off'] keep reference channel in output (only for several
%                  references).
%
% Outputs:
%   dataout     - Input data converted to average reference
%   Chanlocs    - Updated location structure
%   Wout        - ICA weight matrix converted to average reference
%   Sout        - ICA sphere matrix converted to eye()
%   meandata    - (1,dataframes) mean removed from each data frame (point)
%
% Notes: 1) The average reference calculation implements two methods 
%           (from www.egi.com/Technotes/AverageReference.pdf)
%            V'i= (Vi-Vref) - sum(Vi-Vref)/number_of_electrodes
%        2) 'icaweight' conversion of the weight matrix W to average reference:
%        If ica_act = W*data, then data = inv(W)*ica_act; 
%        If R*data is the average-referenced data, 
%        R*data= (R*inv(W))*ica_act and Wout = inv(R*inv(W));
%        The average-reference ICA maps are the columns of inv(Wout).
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

function [data, Elocs, W, S, meandata] = reref(data, ref, varargin)

if nargin<1
  help reref
  return
end
if nargin < 2
    ref = [];
end;

% check inputs
% ------------
g = finputcheck(varargin, { 'icaweight'  'real'    []          [];
                            'method'     'string'  { 'standard' 'withref' }  'standard';
                            'refstate'   { 'string' 'integer' } { { 'common' 'averef' } [1 size(data,1)] }     'common';
                            'exclude'    'integer' [1 size(data,1)]          [];
                            'refloc'     'cell'    []                        {};
                            'keepref'    'string'  {'on' 'off' }             'off';
                            'elocs'      {'integer' 'struct'}  []            [] });
if isstr(g), error(g); end;

chans = size(data,1);
if chans < 2 
  help reref
  return
end
if ref > chans
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
            error('Cannot include old reference if multiple references had been used');
        end;
        disp('Note: old reference electrode included in average reference computation');
        avematrix          = eye(nbchans)-ones(nbchans)*1/(nbchans+1); % reference is a relrevant channel i.e. Cz
        avematrix(end+1,:) = -1/(nbchans+1);                           % potential for the new electrode (previously ref)
        if ~( length(g.elocs) > chans )
            if ~isempty(g.refloc) & ~isempty(g.refloc{1})
                g.elocs(end+1).labels = g.refloc{1};
                g.elocs(end  ).theta  = g.refloc{2};
                g.elocs(end  ).radius = g.refloc{3};
            else
                error('Need a old reference location to include it as a data channel');
            end;
        end;
        rerefchansout(end+1) = chans+1;
    else
        disp('Note: old reference electrode (unless present as data channel) not included in average reference computation');
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
chans = size(avematrix,1);
refmatrix = eye(chans);
if ~isempty(ref)
    
    fprintf('Re-referencing data\n');
    for index = 1:length(ref)
        refmatrix(:,ref(index)) = refmatrix(:,ref(index))-1/length(ref);
    end;
    fprintf('Warning: reference channels have been removed\n');
    if length(ref) > 1 
        if strcmpi(g.keepref, 'off')
            refmatrix([ref g.exclude],:) = []; % supress references and non EEG channels
            refmatrix(:,g.exclude      ) = [];              % supress non EEG channels
            rerefchansout = setdiff(rerefchansout, ref);

            if ~isempty(g.elocs)
                g.elocs(ref) = [];
            end;
        else
            refmatrix(g.exclude,:) = []; % supress non EEG channels
            refmatrix(:,g.exclude) = []; % supress non EEG channels
        end;
    else
        refmatrix([ref g.exclude],:) = []; % supress references and non EEG channels
        refmatrix(:,g.exclude      ) = [];              % supress non EEG channels
        rerefchansout = setdiff(rerefchansout, ref);

        % copy channel location for single ref
        % ------------------------------------
        if ~isempty(g.elocs) & length(g.elocs) > chans & length(ref) == 1
            tmpelocs       = g.elocs(ref);
            g.elocs(ref)   = [];
            g.elocs(end+1) = tmpelocs;
        end;
    end;
end;

% there are faster methods but this one is the simpliest

data(rerefchansout,:) = (refmatrix*avematrix)*data(rerefchans,:); % implement as a matrix multiply
if strcmpi(g.keepref, 'off')
    data = data(setdiff(1:size(data,1), ref),:);
end;

Elocs = g.elocs;
data = reshape(data, size(data,1), dim2, dim3);

% treat optional ica parameters
% -----------------------------
if ~isempty(g.icaweight) 
	winv = pinv(g.icaweight);
    try, 
        W = pinv(refmatrix*avematrix*winv);
	catch,
        error([ 'Weight matrix size is different from the data size, re-referencing impossible' 10 ...
                      '(you have to use the same number of channel in rereferenging (minus excluded ones)' 10 ...
                      'as in the ICA decomposition; best solution is to suppress ICA weights, rereference, then rerun ICA)']);
    end;
    S = eye(size(W,2));
else
    W = []; S = [];
end;
