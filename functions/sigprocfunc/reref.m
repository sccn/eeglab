% reref() - convert common-reference EEG data to average reference
%
% Usage:
%   >> data = reref(data);
%   >> [dataout Elocs Wout Sout meandata] = reref(data, ref, 'key', 'val');
%
% Inputs:
%   data - 2D data matrix (chans,frames*epochs) 
%   ref  - reference channel number: 2 possibilities
%          1) [] - compute average reference. If the 'withref' mode
%             is set, the function recomputes the common reference potential 
%             while averaging.
%          2) 1<= X <= size(data,1): re-reference to channel X. If the 
%             'withref' mode is set, the function compute the potential of
%             of the old common reference.
% 
% Optional inputs:
%   'elocs'     - Electrode location structure.
%   'icaweight' - ICA weight matrix (or weight * sphere)
%   'method'    - ['standard'|'withref'] include reference in the output
%                 data.
%   'refloc'    - reference channel location
%                 Cell array with name and polar coordinates of the channel
%                 { 'label' theta radius }. For 3-D location, include 
%                 the reference as the last channel in the 'elocs' structure.
%
% Outputs:
%   dataout  - Input data converted to average reference
%   Elocs    - Updated location structure
%   Wout     - ICA weight matrix converted to average reference
%   Sout     - ICA sphere matrix converted to eye()
%   meandata - (1,dataframes) mean removed from each data frame (point)
%
% Notes: 1) Average reference calculation implement two methods 
%        from www.egi.com/Technotes/AverageReference.pdf
%        standard: V'i = (Vi - Vref) - sum(Vi-Vref)/number_of_electrode
%        withref:  V'i = (Vi - Vref) - sum(Vi-Vref)/(number_of_electrode+1)
%        The last method also allow to compute the potential of the common
%        reference. Note that the common reference can also be included
%        when re-referencing data.
%        2) If 'icaweight' converts the weight matrix W to average reference:
%        If ica_act = W*data, then data = inv(W)*ica_act; 
%        If R*data is the average-referenced data, 
%        R*data = (R*inv(W))*ica_act and Wout = inv(R*inv(W));
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
if naring < 2
    ref = [];
end;

% check inputs
% ------------
g = finputcheck(varargin, { 'icaweight'  'real'    []          [];
                            'method'     'string'  { 'standard' 'withref' }  [];
                            'refloc'     'cell'    []          [];
                            'elocs'      'struct'  []          [] });
if isstr(g), error(g); end;

chans = size(data,1);
if chans < 2 
  help reref
  return
end

[dim1 dim2 dim3] = size(data);
data = reshape(data, dim1, dim2*dim3);

% compute potential of reference
% ------------------------------
if isempty(ref)
    if strcmp(g.method, 'withref')
        avematrix = eye(chans)-ones(chans)*1/(chans+1);
        avematrix(end+1,:) = -1; % common reference channel
        if ~isempty(g.elocs)
            if (length(g.elocs) == chans) & ~isempty(g.refloc)
                g.elocs(end+1).labels = g.refloc{1};
                g.elocs(end  ).theta  = g.refloc{2};
                g.elocs(end  ).radius = g.refloc{3};
            else 
                error('No location for old common reference, can not introduce it as a new channel');
            end;
        end;        
        meandata = sum(data)/(chans+1);
    else 
        avematrix = eye(chans)-ones(chans)*1/chans;    
        % electrode structure remain unchanged except if reference electrode is present
        if ~isempty(g.elocs) & length(g.elocs) > chans
            g.elocs(end) = [];
        end;
        meandata = sum(data)/chans;
    end;
else % common re-reference
     % -------------------
    if strcmp(g.method, 'withref')
        avematrix = eye(chans);
        avematrix(:,ref) = -1;
        if ~isempty(g.elocs)
            if length(g.elocs) > chans
                tmpelocs     = g.elocs(ref);
                g.elocs(ref) = g.elocs(end);
            elseif ~isempty(g.refloc)
                g.elocs(end+1) = g.elocs(ref);
                g.elocs(ref).labels = g.refloc{1};
                g.elocs(ref).theta  = g.refloc{2};
                g.elocs(ref).radius = g.refloc{3};
            else 
                error('No location for old common reference, can not introduce it as a new channel');
            end;
        end;
    else
        avematrix = eye(chans);
        avematrix(:,ref) = -1;
        avematrix(ref,:) = [];
        if ~isempty(g.elocs)
            tmpelocs     = g.elocs(ref);
            g.elocs(ref:end-1) = g.elocs(ref+1:end);
            g.elocs(end)       = tmpelocs;
        end;
    end;
end;
    
data = avematrix*data; % implement as a matrix multiply
% there are faster methods but this one is the simpliest
Elocs = g.elocs;

% treat optional ica parameters
% -----------------------------
if ~isempty(g.icaweight) 
	winv = pinv(g.icaweigth);
    try, 
        W = pinv(avematrix*winv);
	catch,
        error('Weight matrix size is different from the data size, re-referencing impossible');
    end;
    S = eye(chans);
end;
