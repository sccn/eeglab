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
%   'method'     - ['standard'|'withref'] Include reference channel in output. 
%   'refstate  ' - ['common'|'averef'] Current average reference state.
%                  Use this parameter to re-reference data to a given channel if data
%                  is already average referenced (by setting it to 'averef') 
%                  Default is 'common'.
%   'refloc'     - Reference channel location -- a cell array with name and 
%                  polar coordinates of the channel, { 'label' theta radius }. 
%                  For 3-D location, include the reference as the last channel 
%                  in the 'chanlocs' structure.
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
%            V'i = (Vi-Vref) - sum(Vi-Vref)/number_of_electrodes
%        2) 'icaweight' conversion of the weight matrix W to average reference:
%        If ica_act = W*data, then data = inv(W)*ica_act; 
%        If R*data is the average-referenced data, 
%        R*data=(R*inv(W))*ica_act and Wout = inv(R*inv(W));
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
                            'method'     'string'  { 'standard' 'withref' }  '';
                            'refstate'   'string'  { 'common' 'averef' }   'common';
                            'refloc'     'cell'    []          {};
                            'elocs'      {'integer' 'struct'}  []          [] });
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

% compute potential of reference
% ------------------------------
if isempty(ref)
    avematrix = eye(chans)-ones(chans)*1/chans;
    if strcmp(g.refstate, 'averef')
        fprintf('Undoing average reference\n');
        avematrix = pinv(avematrix);
    else 
        fprintf('Computing average reference\n');
    end;
    % electrode structure remain unchanged except if reference electrode is present
    if ~isempty(g.elocs) & length(g.elocs) > chans
        g.elocs(end) = [];
    end;
    meandata = sum(data)/chans;
else % common re-reference
     % -------------------
    
    % compute the inverse average transformation matrix
    % -------------------------------------------------
    if strcmp(g.refstate, 'averef')
        fprintf('Undoing average reference\n');
        invavematrix = eye(chans)-ones(chans)/chans;        
    else 
        invavematrix = [];
    end;
    fprintf('Re-referencing data\n');
    if strcmp(g.method, 'withref')
        if length(ref) > 1
            % dealing with multiple references
            % --------------------------------
            avematrix = eye(chans);
            for index = 1:length(ref)
                avematrix(:,ref) = -1/length(ref);
            end;
            fprintf('Warning: reference channels have been removed');
            avemaxtrix(ref(2:end),:) = []; % supress references
            if ~isempty(g.elocs)
                g.elocs(ref(2:end)) = [];
            end;
            ref = ref(1);
            if length(g.elocs) > size(avematrix,1)
                g.elocs(ref) = g.elocs(end);
                g.elocs(end) = [];
            elseif ~isempty(g.refloc) & ~isempty(g.refloc{1})
                g.elocs(ref).labels = g.refloc{1};
                g.elocs(ref).theta  = g.refloc{2};
                g.elocs(ref).radius = g.refloc{3};
            else 
                error('No location for old common reference, cannot introduce it as a new channel');
            end;            
        else
            % dealing with a single ref. channel
            % ----------------------------------
            avematrix = eye(chans);
            avematrix(:,ref) = -1;
            if ~isempty(g.elocs)
                if length(g.elocs) > chans
                    tmpelocs     = g.elocs(ref);
                    g.elocs(ref) = g.elocs(end);
                    g.elocs(end) = tmpelocs;
                elseif ~isempty(g.refloc)
                    g.elocs(end+1) = g.elocs(ref);
                    g.elocs(ref).labels = g.refloc{1};
                    g.elocs(ref).theta  = g.refloc{2};
                    g.elocs(ref).radius = g.refloc{3};
                else 
                    error('No location for old common reference, cannot introduce it as a new channel');
                end;
            end;
        end;
    else
        if length(ref) > 1
            % dealing with multiple references (do not include reference)
            % --------------------------------
            avematrix = eye(chans);
            for index = 1:length(ref)
                avematrix(:,ref) = -1/length(ref);
            end;
            fprintf('Warning: reference channels have been removed');
            avemaxtrix(ref,:) = []; % supress references
            if ~isempty(g.elocs)
                g.elocs(ref) = [];        
            end;    
        else
            % dealing with a single ref. channel
            % ----------------------------------
            avematrix = eye(chans);
            avematrix(:,ref) = -1;
            avematrix(ref,:) = [];
            if ~isempty(g.elocs)
                g.elocs(end)       = g.elocs(ref);
                g.elocs(ref)       = [];
            end;
        end;
    end;
    
    if ~isempty(invavematrix)
        avematrix = avematrix * pinv(invavematrix);
    end;
end;
    
data = avematrix*data; % implement as a matrix multiply
% there are faster methods but this one is the simpliest
Elocs = g.elocs;
data = reshape(data, size(data,1), dim2, dim3);

% treat optional ica parameters
% -----------------------------
if ~isempty(g.icaweight) 
	winv = pinv(g.icaweight);
    try, 
        W = pinv(avematrix*winv);
	catch,
        error('Weight matrix size is different from the data size, re-referencing impossible');
    end;
    S = eye(size(W,2));
end;
