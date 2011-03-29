% pop_fmrib_fastr() -  GUI for fmrib_fastr: Remove gradient
%   artifacts from EEG data collected during FMRI using FMRI Artifact
%   Slice Template Removal [Niazy 2006] and reduce residuals by removing
%   artifact principal components.
%
% This program prompts the user for settings required to carry out
% removal of gradient artifacts from EEG data by running fmrib_fastr.m.
%
% Usage:
%   >> [EEGOUT, COM] = pop_fmrib_fastr(EEG) - 
%
% Command Line usage:
%   >> [EEGOUT, command] = pop_fmrib_fastr(EEG,lpf,L,Win,etype,...
%         strig,anc_chk,trig_correct,Volumes,Slices,pre_frac,exc_chan,NPC); 
%
% Inputs: (see details below and fmrib_fastr.m)
%       EEG - EEG structure.
%       lpf - Low pass filter cutoff (default: [ ]=70).
%       L - Interpolation folds (default: [ ]=10).
%       Win - Number of artifacts in avg. window (default: [ ]=30).
%       etype - Name of FMRI slice (slice timing) event.  Unlike
%           fmrib_fastr.m, this is the name of the event.  fmrib_fastr.m 
%           takes a vector with slice locations.
%       strig:          1 for slice triggers (default)
%                       0 for volume/section triggers
%       anc_chk         1 to perform ANC
%                       0 No ANC
%       trig_correct:   1 to correct for missing triggers;
%                       0 to NOT correct.
%                       (default: [ ]=0)
%       Volumes:        FMRI Volumes.  Needed if trig_correct=1; 
%                       otherwise use [ ];
%       Slices:         FMRI Slices/Vol. usage like Volumes.
%       pre_frac:       Relative location of slice triggers to actual start
%                       of slice (slice artifact).  Value between 0 - 1.
%                       0 = slice trigger (event) points to absolute start
%                       of slice (artifact); 1 = points to absolute end.
%                       (default: [ ]=0.03).
%       exc_chan:       Channels to exclude from residual artifact
%                       principle component fitting (OBS).  Use for EMG, 
%                       EOG or other non-eeg channels.
%       NPC:            Number of principal components to fit to residuals.
%                       0 to skip OBS fitting and subtraction.
%                       'auto' for automatic order selection (default)
%
% GUI Inputs:
%------------
% Low pass filter cutoff frequency (HZ): Low pass filtering is required to
%     carry out FASTR.  Enter the desired cutoff frequency. Default: 0 Hz. 
%     (No Filtering).
%
% Interpolation (up-sampling) folds:  in order to correctly remove the 
%     gradient artifacts, the algorithm up-samples the data.  It is 
%     recommended that you ups-ample the data to bring the sampling 
%     frequency to about 20 kHz. The default up-sampling factor is 10.
%
% Averaging window size: In here enter how many realization of the slice 
%     artifacts are to be included in the averaging window.  The bigger the 
%     window the better the approximation of the true artifact is, but the 
%     slower the adaptation is to any changes in the shape of the artifact 
%     (due to movement for example).  
%
% Artifact-timing event type:  The program requires that the data in EEGLAB 
%    have events that correspond to the timing of each FMRI repeated artifact
%    period.  Select between slice or volume/event triggers.
%    The pop-up menu displays all available event types.  
%    Select the one corresponding to the slice timing.  
%
% Adaptive Noise Cancellation (ANC):  This checkbox turns on ANC.  This
%    will make the entire process much slower.  To use ANC, a low-pass 
%    filter must be specified.  If not, an LPF of 70Hz will be used.
%
% Correct for missing triggers:  Check this option if you want FASTR to add 
%    missed triggers due technical faults.  The algorithm that corrects for 
%    this needs to know the number of FMRI volumes and Slices collected 
%    during the acquisition of the EEG data.  It then calculated the 
%    theoretical number of correct triggers to expect, and add any missing 
%    ones.  This feature has not been tested extensively, so use it with 
%    caution.  The draw back of not using it is that you could have 
%    segments of artifact-contaminated data where the slice events are 
%    missing.
%  
% Non-EEG Channels:  Enter any channels in the data that are not EEG 
%    channels.  This is because RAPCO (see fmrib_fastr.m) might remove
%    features such as QRS complexes in ECG data.  The artifacts will still
%    removed but the residuals will not be further processed with RAPCO.
%
% Advanced Options: 
%     
%     Relative trigger location: During acquisition of each FMRI slice, the 
%     MRI machine fires a trigger.  The time at which the trigger is fired 
%     relative to the start of the slice acquisition may vary from one MRI 
%     machine to another.  In this area, you can define the position of the 
%     trigger relative to the start of the slice.  A value of zero 
%     indicates the trigger is at the exact beginning of the slice 
%     collection and 1 indicates at the end.  Usually the trigger is much 
%     closer to the beginning than the end.  The default value is 0.03. 
%
%     Number of residuals principal components to remove:FASTR (by default 
%     - 'auto') automatically determines the number of residual PCs to 
%     remove using the plot of ordered eigenvalus.  Otherwise, one can 
%      simply enter the number of residual PCs to remove. Tor different 
%     data sets originally sampled at different frequencies, different 
%     number of PCs will be needed. The best way to know is to do a PCA 
%     on a martix of the slice artifacts (after subtracting the artifact 
%     but before filtering) and look at the eigenvalue plot.  If you are 
%     not sure what that means or how to do this use the default or 
%     experiment with different number of PCs starting with 1 and 
%     increasing with steps of 1. NOTE: TO SKIP OBS FITTING AND 
%     SUBTRACTION, SET THIS = 0.
%
%
% Author: Rami Niazy, FMRIB Centre, University of Oxford.
%
% 
% Copyright (c) 2004 University of Oxford

% Copyright (C) 2004 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
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

% SEP 16, 2005
% Help update

% AUG 5, 2005
% ANC made optional.
% Handles 'auto' order selection of PCs.

% 12 MAR 2005
% Program edited to work with volume/section triggers.

% JAN 14, 2005
% Updated help.
% Fixed history command output

% DEC 23, 2004
% update (c)

% updated: Dec 14, 2004.

function [EEGOUT, command] = pop_fmrib_fastr(EEG,lpf,L,Win,etype,...
    strig,anc_chk,trig_correct,Volumes,Slices,pre_frac,exc_chan,NPC); 

EEGOUT = [];
command = '';

% check inputs
%--------------------

if nargin < 1
    error('Not enough input arguments','fmrib_fastr() error!');
elseif nargin > 1
    if nargin < 11
        error('Not enough input arguments','fmrib_fastr() error!');
    end
end

if isempty(EEG.data)
    errordlg2(strvcat('Error: no data'), 'Error');
    error('pop_fmrib_fastr error: no data','fmrib_fastr() error!'); return;
end

% check for DSP Toolbox
%----------------------
if ~exist('firls')
   error('FASTR requires the signal processing toolbox.','fmrib_fastr() error!');
end

%Launch user input GUI and read settings
%--------------------------------------
if nargin==1
    tmpevent = EEG.event;
    etypes=unique({tmpevent.type});
    guifig=guifastr(etypes);
    p=guidata(guifig);
else
    p.opstatus=1;
    
    p.param.etype=etype;
    
    if isempty(lpf)
        p.param.lpf='0';
    else
        p.param.lpf=num2str(lpf);
    end
    
    if isempty(strig)
        p.param.strig=1;
    else
        p.param.strig=strig;
    end
    
    if isempty(anc_chk)
        p.param.anc_chk=0;
    else
        p.param.anc_chk=anc_chk;
    end
    
    if isempty(L)
        p.param.L='10';
    else
        p.param.L=num2str(L);
    end
    
    if isempty(Win)
        p.param.window='30';
    else
        p.param.window=num2str(Win);
    end
    
    if isempty(trig_correct) | trig_correct==0
        p.param.trigcorrect_check=0;
        p.param.volumes='0';
        p.param.slices='0';
    elseif trig_correct==1
        p.param.trigcorrect_check=1;
        if Volumes<=0 | Slices<=0 | isempty(Volumes) | isempty(Slices)
            error('Invalid Slices and Volumes entries for trigger correction',...
                'fmrib_fastr() error!');
        end
        p.param.volumes=num2str(Volumes);
        p.param.slices=num2str(Slices);
    else
        error('trig_correct input should be 1 or 0.  See help.',...
            'fmrib_fastr() error!');
    end    
    
    if isempty(pre_frac)
        p.param.pre_frac='0.03';
    else
        p.param.pre_frac=num2str(pre_frac);
    end
    
    if isempty(exc_chan)
        p.param.exc='';
    else
        p.param.exc=num2str(exc_chan);
    end
    
    if isempty(NPC)
        p.param.npc='auto';
    else
        if ischar(NPC)
            p.param.npc='auto';
        else
            p.param.npc=num2str(NPC);
        end
    end
end
        
    

%process using settings
%----------------------------------
switch p.opstatus
case 0 % this is cancel button setting
    delete(guifig);
case 1 % this is ok button setting
    if nargin==1
        delete(guifig);
        etype=char(etypes{p.param.etype});
    end
    
    Trigs=[];
	for E=1:length(EEG.event)
        if strcmp(EEG.event(E).type,etype)
            Trigs(end+1)=round(EEG.event(E).latency);
        end
    end
    
    comline = sprintf('%s',p.param.lpf);
	comline = sprintf('%s,%s,%s,Trigs',comline,p.param.L,p.param.window);
    comline = sprintf('%s,%s,%s',comline,num2str(p.param.strig),num2str(p.param.anc_chk));
    comline = sprintf('%s,%s,%s,%s', comline, num2str(p.param.trigcorrect_check),p.param.volumes,p.param.slices);
    if strcmp(lower(p.param.npc),'auto')
        comline = sprintf('%s,%s,[%s],''%s''', comline,p.param.pre_frac,p.param.exc,p.param.npc);
    else
        comline = sprintf('%s,%s,[%s],%s', comline,p.param.pre_frac,p.param.exc,p.param.npc);
    end
    comline = [ 'fmrib_fastr(EEG,' comline ');'];
    EEGOUT=eval(comline);
    fprintf('-----------------------------------\n');
    
    command = sprintf('%s',p.param.lpf);
    command = sprintf('%s,%s,%s,''%s''',command,p.param.L,p.param.window,etype);
    command = sprintf('%s,%s,%s',command,num2str(p.param.strig),num2str(p.param.anc_chk));
    command = sprintf('%s,%s,%s,%s', command, num2str(p.param.trigcorrect_check),p.param.volumes,p.param.slices);
    if strcmp(lower(p.param.npc),'auto')
        command = sprintf('%s,%s,[%s],''%s''', command,p.param.pre_frac,p.param.exc,p.param.npc);
    else
        command = sprintf('%s,%s,[%s],%s', command,p.param.pre_frac,p.param.exc,p.param.npc);
    end
    command = [ 'EEG = pop_fmrib_fastr(EEG,' command ');'];
end
return;
