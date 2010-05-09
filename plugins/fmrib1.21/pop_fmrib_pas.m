% pop_fmirb_pas() -  GUI for fmrib_pas: Remove hear-related pulse artifacts
% from EEG data.
%
% This program prompts the user for settings required to carry out
% removal of pulse artifacts from EEG data by running fmrib_pas.m.
%
% Usage:
%    >> [EEGOUT, COM] = pop_fmrib_pas(EEG)    %  this launches gui
% or >> [EEGOUT, COM] = pop_fmrib_pas(EEG,qrseventtype,method) 
% or >> [EEGOUT, COM] = pop_fmrib_pas(EEG,qrseventtype,method,npc)  for
% command  line usage.
%
% Inputs:
% qrseventtype: the name of the heart-beat/QRS events in EEG
% method: 'obs', 'mean' , 'gmean', or 'median'. These correspond to the 
% three methods below
%
% GUI Inputs:
%------------
% QRS / Heart Beat Events:  This popup menu displays available event types
%   in the data.  Select the event that corresponds to QRS complex locations
%   from ECG data.  Also see pop_fmrib_qrsdetect.
%
% Artifact Template Formation Method:
%   *Optimal Basis Set: Does a PCA on a martix of all
%   the heart artifacts then fits the first N components to each artifact.
%   The default number of components is 4.  
%   *Simple Mean:  Simply averages successive pulse artifacts
%   *Gaussian-Weighted Mean:  Averages artifacts after multiplying by
%       a Gaussian window weights to emphasise current artifact shape
%       and reduce effect of further artifacts.
%   *Median:  Uses a median filter of artifacts to form template.
%  
%
% Author: Rami Niazy, FMRIB Centre, University of Oxford.
%
% 
%
% Copyright (c) 2004 University of Oxford.
%


% Copyright (C) 2004 University of Oxford
% Author:  Rami K. Niazy, FMRIB Centre
%          rami@fmrib.ox.ac.uk
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

% JAN 10, 2004
% Fixed Command History Output

% DEC 23, 2004
% Update (c)

% Dec 15, 2004
% Handles 'obs' method and
% passes N of PCs to fmrib_pas.

% Oct 26, 2004
% Function now processes events
% and passes vector of peak locations
% to fmrib_pas

function [EEGOUT, command] = pop_fmrib_pas(EEG,varargin); 

EEGOUT = [];
command = '';

% check inputs
%--------------------
if nargin > 1
    if nargin < 3
        error('Incorrect number of input arguments');
    end
elseif nargin < 1
    error('Incorrect number of input arguments');
end

if isempty(EEG.data)
    errordlg2(strvcat('Error: no data'), 'Error');
    error('pop_fmrib_fastr error: no data','fmrib_fastr() error!'); return;
end

% check for DSP Toolbox
%----------------------
if ~exist('firls')
   error('fmrib_pas requires the signal processing toolbox.','fmrib_pas() error!');
end

%Launch user input GUI for 1 input or assign user inputs
%--------------------------------------------------------
if nargin==1
    etypes=unique({EEG.event.type});
	guifig=guipas(etypes);
else
    p.opstatus=1;
    etype=varargin{1};
    p.param.method=varargin{2};
    p.param.npc='4';
    switch p.param.method
    case 'obs'
        if nargin==4
            p.param.npc=num2str(varargin{3});
        end
    case 'mean'
    case 'gmean'
    case 'median'
    otherwise 
        error('Unknown method','pop_fmrib_pas()! error.');
    end        
end

%when function returns read settings
%----------------------------------
if nargin==1
    p=guidata(guifig);
end

switch p.opstatus
case 0 % this is cancel button setting
    delete(guifig);
case 1 % this is ok button setting
    
    if nargin==1
        delete(guifig);
        etype=char(etypes{p.param.qrsevent});
    end
    QRSevents=[];
	for E=1:length(EEG.event)
        if strcmp(EEG.event(E).type,etype)
            QRSevents(end+1)=round(EEG.event(E).latency);
        end
	end
    
    comline = sprintf('QRSevents,''%s'',%s',p.param.method,p.param.npc);
    comline = [ 'fmrib_pas(EEG,' comline ');'];
    EEGOUT=eval(comline);
    command = sprintf('''%s'',''%s'',%s',etype,p.param.method,p.param.npc);
    command=['pop_fmrib_pas(EEG,' command ');'];
end
command = sprintf('EEG = %s',command);
return;

