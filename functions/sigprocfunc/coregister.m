% coregister() -  Co-register measured or template electrode locations with a 
%                 a reference channel locations file. For instance if you
%                 want to perform dipole modeling you have to coregister
%                 (align) your channel electrodes with the model (and the
%                 easiest way to do that is to coregister your channel
%                 electrodes with the electrodes file associated with the
%                 model. To use coregister(), one may for instance use the default 
%                 MNI head and 10-5 System locations from Robert Oostenveld, used in  
%                 the dipfit2() dipole modeling function as a reference. 
%                 Use coregister() to linearly align or nonlinearly warp 
%                 subsequently recorded or nominally identified (e.g., 'Cz') sets of 
%                 channel head locations to the reference locations. 
%                 Both channel locations and/or fiducial (non-electrode) 
%                 locations can then be used by coregister() to linearly align 
%                 or nonlinearly warp a given or measured montage to the reference 
%                 locations. In its (default) manual mode, coregister() produces 
%                 an interactive gui showing the imported and reference channel 
%                 locations on the imported head mesh (if any), allowing the user 
%                 to make additional manual adjustments using gui text entry boxes, 
%                 and to rotate the display using the mouse.
% Usage:
%        >>  coregister(EEG.chanlocs); % Show channel locations in the coregister() 
%                 % gui, for align, warp, or manual mode co-registration with the 
%                 % dipfit2() model head mesh and 10-5 System reference locations.
%                 % Note: May need to scale channel x,y,z positions by 100 to see 
%                 % the imported locations in the default dipfit2() head mesh.
%
%        >> [chanlocs_out transform] = coregister( chanlocs, reflocs, 'key', 'val' )
%                 % Perform custom co-registration to a reference locations and
%                 % (optional) head mesh. If a ('warp' or 'alignfid') mode is
%                 % specified, no gui window is produced.
% Inputs:
%    chanlocs   - (EEG.chanlocs) channel locations structure (or file) to align
%                 For file structure, see >> headplot example 
%                                      or >> headplot cartesian
%    reflocs    - reference channel locations structure (or file) associated 
%                 with the head mesh ('mesh' below) {default|[] -> 10-5 System locs}
%
% Optional 'keywords' and arguments ([ ]):
%    'mesh'     - [cell array|string] head mesh. May contain {points triangles} 
%                 or {points triangles normals} (see >> help plotmesh 
%                 for details). May also be the name of a head mesh .mat
%                 file (several formats recognized). By default, loads the 
%                 dipfit2() MNI head mesh file. Compare 'mheadnew.mat', the 
%                 mesh used by headplot(); 'mesh',[] shows no head mesh.
%   'warpmethod' - ['rigidbody'|'globalrescale'|'traditional'|'nonlin1'|
%                 'nonlin2'|'nonlin3'|'nonlin4'|'nonlin5']
%                 'traditional' calls the dipfit2.* function traditionaldipfit()
%                 all others are enacted by electrodenormalize()
%                 {default: 'traditional}
%   'transform' - [real array] homogenous transformation matrix (>>help 
%                 electrodenormalize) or a 1x9 matrix containing traditional 
%                 9-parameter "Talairach model" transformation (>> help traditional) 
%                 used to calculate locs_out.
%   'chaninfo1' - [EEG.chaninfo struct] channel information structure for the 
%                 montage to be coregistered. May contain (no-electrode) fiducial
%                 locations.
%   'chaninfo2' - [EEG.chaninfo struct] channel information structure for 
%                 the reference montage.
%   'helpmsg'   - ['on'|'off'] pop-up help message when calling function.
%                 {default: 'off'}
%   'showlabels1' - ['on'|'off'] show channel labels for first montage. 
%                 Default is 'off'.
%   'showlabels2' - ['on'|'off'] show channel labels for second montage. 
%                 Default is 'off'.
%
% Optional 'keywords' for MANUAL MODE (default):
%   'manual'    - ['on'|'off'|'show'] Pops up the coregister() gui window to 
%                 allow viewing the current alignment, performing 'alignfid' or 
%                 'warp' mode co-registration,  and making manual
%                 adjustments. Default if 'on'. 'off' does not pop up any 
%                 window and 'show' does not allow to change coregistration.
%
% Optional 'keywords' for ALIGN MODE:
%    'alignfid' - {cell array of strings} = labels of (fiducial) channels to use 
%                 as common landmarks for aligning the new channel locations to 
%                 the reference locations. These labels must be in both channel 
%                 locations (EEG.chanlocs) and/or EEG.chaninfo structures (see 
%                 'chaninfo1' and 'chaninfo2' below). Will then apply a linear 
%                 transform including translation, rotation, and/or uniform 
%                 scaling so as to best match the two sets of landmarks. 
%                 See 'autoscale' below.
%   'autoscale' - ['on'|'off'] autoscale electrode radius when aligning 
%                 fiducials. {default: 'on'}
%
% Optional 'keywords' for WARP MODE:
%    'warp'     - {cell array of strings} channel labels used for non-linear 
%                 warping of the new channel locations to reference locations.
%                 For example: Specifying all 10-20 channel labels will
%                 nonlinearly warp the 10-20 channel locations in the new
%                 chanlocs to the 10-20 locations in the reference chanlocs.
%                 Locations of other channels in the new chanlocs will be 
%                 adjusted so as to produce a smooth warp. Entering the 
%                 string 'auto' will automatically select the common channel 
%                 labels.
%
% Outputs:
%  chanlocs_out - transformed input channel locations (chanlocs) structure
%  transform    - transformation matrix. Use traditionaldipfit() to convert
%                 this to a homogenous transformation matrix used in
%                 3-D plotting functions such as headplot().
%
% Note on how to create a template:
%                 (1) Extract a head mesh from a subject MRI (possibly in 
%                 Matlab using the function isosurface). 
%                 (2) Measure reference locations on the subject's head. 
%                 (3) Align these locations to the extracted subject head mesh 
%                 (using the coregister() graphic interface). The aligned locations 
%                 then become reference locations for this mesh.
%
% Example:
%   % This return the coregistration matrix to the MNI brain
%   dipfitdefs;
%   [newlocs transform] = coregister(EEG.chanlocs, template_models(2).chanfile, ...
%        'warp', 'auto', 'manual', 'off');
%
%   % This pops-up the coregistration window for the BEM (MNI) model
%   dipfitdefs;
%   [newlocs transform] = coregister(EEG.chanlocs, template_models(2).chanfile, ...
%        'mesh', template_models(2).hdmfile);
%
% See also: traditionaldipfit(), headplot(), plotmesh(), electrodenormalize(). 
%
% Note: Calls Robert Oostenveld's FieldTrip coregistration functions for
%       automatic coregistration.
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2005-06

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2005, arno@salk.edu
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

% Here we bypass channel info
%    'chaninfo'  - [struct] channel information structure for the first
%                  dataset (may use a cell array { struct struct } to enter
%                  channel location info for both channel location struct.
%

function [ chanlocs1, transformmat ] = coregister(chanlocs1, chanlocs2, varargin)

if nargin < 1
  help coregister
  error('Not enough inputs');
end

if nargin < 2
    chanlocs2 = [];
end
    % manually for s1:  transf = [4 0 -50 -0.3 0 -1.53 1.05 1.1 1.1]
    % manually for s2:  transf = [-4 -6 -50 -0.37 0 -1.35 1.1 1.15 1.1]

% if isempty(chanlocs2)
%     dipfitdefs; chanlocs2 = template_models(1).chanfile; 
% end

% undocumented commands run from GUI
% ----------------------------------
if ischar(chanlocs1) 
    if ~strcmpi(chanlocs1, 'redraw') && ~strcmpi(chanlocs1, 'fiducials') && ~strcmpi(chanlocs1, 'warp')
        chanlocs1 = readlocs(chanlocs1);
    else
        com = chanlocs1;
        fid = chanlocs2;

        % update GUI
        % ----------
        if strcmpi(com, 'redraw'), redrawgui(fid); return; end

        % select electrodes and warp montage
        % ----------------------------------
        dat = get(fid, 'userdata');
        if strcmpi(com, 'fiducials')
            [clist1, clist2] = pop_chancoresp( dat.elec1, dat.elec2, 'autoselect', 'fiducials');
            try
                [ ~, transform ] = align_fiducials(dat.elec1, dat.elec2, dat.elec1.label(clist1), dat.elec2.label(clist2));
                if ~isempty(transform), dat.transform = transform; end
            catch
                warndlg2(strvcat('Transformation failed', lasterr));
            end
        elseif strcmpi(com, 'warp')
            [clist1, clist2] = pop_chancoresp( dat.elec1, dat.elec2, 'autoselect', 'all');

            % copy electrode names
            if ~isempty(clist1)
                tmpelec2 = dat.elec2;
                for index = 1:length(clist2)
                    tmpelec2.label{clist2(index)} = dat.elec1.label{clist1(index)};
                end
                %try,
                    [ ~, dat.transform ] = warp_chans(dat.elec1, tmpelec2, tmpelec2.label(clist2), 'traditional');
                %catch,
                %    warndlg2(strvcat('Transformation failed', lasterr));
                %end
            end
        end
        set(fid, 'userdata', dat);
        redrawgui(fid); 
        return;
    end
end

% check input arguments
% ---------------------
% defaultmesh = 'D:\matlab\eeglab\plugins\dipfit2.0\standard_BEM\standard_vol.mat';
defaultmesh = 'standard_vol.mat';
g = finputcheck(varargin, { 'alignfid'   'cell'  {}      {};
                            'warp'       { 'string','cell' }  { {} {} }      {};
                            'warpmethod' 'string'  {'rigidbody', 'globalrescale', 'traditional', 'nonlin1', 'nonlin2', 'nonlin3', 'nonlin4', 'nonlin5'} 'traditional';
                            'chaninfo1'  'struct' {}    struct('no', {}); % default empty structure
                            'chaninfo2'  'struct' {}     struct('no', {}); 
                            'transform'  'real'   []      [];
                            'manual'     'string' { 'on','off','show' } 'on'; % -> pop up window
                            'showlabels1'  'string' { 'on','off' } 'off';
                            'showlabels2'  'string' { 'on','off' } 'off';
                            'title'        'string' { } '';
                            'autoscale'  'string' { 'on','off' } 'on';
                            'helpmsg'    'string' { 'on','off' } 'off';
                            'mesh'       ''      []   defaultmesh });
if ischar(g), error(g); end

% load mesh if any
% ----------------
if ~isempty(g.mesh)
    if ischar(g.mesh)
        try
            g.orimesh = g.mesh;
            g.mesh  = load(g.mesh);
        catch, g.mesh = [];
        end
    end
    if ischar(g.mesh)
      if ~exist(g.mesh,'file') 
        fprintf('coregister(): mesh file not found\n');
      end
    end
    if ~isempty(g.mesh)
        if isstruct(g.mesh)
            if isfield(g.mesh, 'SurfaceFile') % Brainstrom leadfield
                p = fileparts(fileparts(fileparts(g.orimesh)));
                try
                    g.mesh = load('-mat', fullfile(p, 'anat', g.mesh.SurfaceFile));
                catch
                    error('Cannot find Brainstorm mesh file')
                end
            end
                
            if isfield(g.mesh, 'vol')
                if isfield(g.mesh.vol, 'r')
                    [X, Y, Z] = sphere(50);
                    dat.meshpnt = { X*max(g.mesh.vol.r) Y*max(g.mesh.vol.r) Z*max(g.mesh.vol.r) };
                    dat.meshtri = [];
                else
                    dat.meshpnt = g.mesh.vol.bnd(1).pnt;
                    dat.meshtri = g.mesh.vol.bnd(1).tri;
                end
            elseif isfield(g.mesh, 'bnd')
                dat.meshpnt = g.mesh.bnd(1).pnt;
                dat.meshtri = g.mesh.bnd(1).tri;
            elseif isfield(g.mesh, 'TRI1')
                dat.meshpnt = g.mesh.POS;
                dat.meshtri = g.mesh.TRI1;
            elseif isfield(g.mesh, 'vertices')
                dat.meshpnt = g.mesh.vertices;
                dat.meshtri = g.mesh.faces;
            elseif isfield(g.mesh, 'Vertices')
                dat.meshpnt = g.mesh.Vertices;
                dat.meshtri = g.mesh.Faces;
            else
                error('Unknown Matlab mesh file');
            end
        else
            dat.meshpnt = g.mesh{1};
            dat.meshtri = g.mesh{2};
        end
    else
        dat.meshpnt = [];
        dat.meshtri = [];
    end
else
    dat.meshpnt = [];
    dat.meshtri = [];
end

% transform to arrays chanlocs1
% -------------------------
TMP          = eeg_emptyset;
TMP.chanlocs = readlocs(chanlocs1);
TMP.chaninfo = g.chaninfo1;
TMP.nbchan = length(TMP.chanlocs);
cfg   = eeglab2fieldtrip(TMP, 'chanloc_withfid');
elec1 = cfg.elec;
if isfield(elec1, 'elecpos')
    elec1.pnt = elec1.elecpos;
end

% transform to arrays chanlocs2
% -------------------------
if ~isempty(chanlocs2)
    TMP   = eeg_emptyset;
    [TMP.chanlocs tmp2 tmp3 ind1] = readlocs(chanlocs2);
    TMP.chaninfo                  = g.chaninfo2;
    TMP.nbchan = length(TMP.chanlocs);
    cfg   = eeglab2fieldtrip(TMP, 'chanloc_withfid');
    elec2 = cfg.elec;
else 
    elec2 = [];
    dat.transform = [ 0 0 0 0 0 0 1 1 1 ];
end

% copy or compute alignment matrix
% --------------------------------
if ~isempty(g.transform)
    dat.transform = g.transform;
elseif ~isempty(elec2)

    % perfrom alignment
    % -----------------
    if strcmpi(g.autoscale, 'on')
        avgrad1 = sqrt(sum(elec1.pnt.^2,2));
        avgrad2 = sqrt(sum(elec2.pnt.^2,2));
        ratio   = mean(avgrad2)/mean(avgrad1);
    else
        ratio = 1;
    end
    if ~isempty(g.alignfid)
        
        % autoscale
        % ---------
        [ electransf transform ] = align_fiducials(electmp, elec2, g.alignfid);
        if ~isempty(transform), dat.transform = [ transform(1:6)' ratio ratio ratio ]; end
        
    elseif ~isempty(g.warp)
        if ischar(g.warp)
            [clist1 clist2] = pop_chancoresp( elec1, elec2, 'autoselect', 'all', 'gui', 'off');
            % copy electrode names
            if isempty(clist1)
                disp('Warning: cannot wrap electrodes (no common channel labels)');
            else
                tmpelec2 = elec2;
                for index = 1:length(clist2)
                    tmpelec2.label{clist2(index)} = elec1.label{clist1(index)};
                end
                try
                    [ ~, dat.transform ] = warp_chans(elec1, tmpelec2, tmpelec2.label(clist2), 'traditional');
                catch
                    warndlg2(strvcat('Transformation failed', lasterr));
                end
            end
        else
            [ ~, dat.transform ] = warp_chans(elec1, elec2, g.warp, g.warpmethod);
        end
    else
        dat.transform = [0 0 0 0 0 0 ratio ratio ratio];
    end
    
end

% manual mode off
% ---------------
if strcmpi(g.manual, 'off') 
    transformmat = dat.transform;
    dat.elec1    = elec1;
    if size(dat.transform,1) > 1
        dat.electransf.pnt = dat.transform*[ dat.elec1.pnt ones(size(dat.elec1.pnt,1),1) ]';
    else
        dat.electransf.pnt = traditionaldipfit(dat.transform)*[ dat.elec1.pnt ones(size(dat.elec1.pnt,1),1) ]';
    end
    dat.electransf.pnt   = dat.electransf.pnt(1:3,:)';
    dat.electransf.label = dat.elec1.label;
    chanlocs1    = dat.electransf;
    return; 
end

% find common electrode names
% ---------------------------
dat.elec1      = elec1;
dat.elec2      = elec2;
dat.elecshow1  = 1:length(elec1.label);
if ~isempty(elec2)
    dat.elecshow2  = 1:length(elec2.label);
else
    dat.elecshow2  = [];
end
dat.color1     = [0 1 0];
dat.color2     = [1 .75 .65]*.8;
%dat.color2     = [1 0 0];
dat.label1     = strcmpi(g.showlabels1, 'on');
dat.label2     = strcmpi(g.showlabels2, 'on');
dat.meshon     = 1;
fid = figure('userdata', dat, 'name', 'coregister()', 'numbertitle', 'off');
try, icadefs; catch, end

if strcmpi(g.manual, 'on') 
    header    = 'dattmp = get(gcbf, ''userdata'');';
    footer    = 'set(gcbf, ''userdata'', dattmp); clear dattmp; coregister(''redraw'', gcbf);';
    cbright   = [ header 'dattmp.transform(1) = str2num(get(gcbo, ''string''));' footer ];
    cbforward = [ header 'dattmp.transform(2) = str2num(get(gcbo, ''string''));' footer ];
    cbup      = [ header 'dattmp.transform(3) = str2num(get(gcbo, ''string''));' footer ];
    cbpitch   = [ header 'dattmp.transform(4) = str2num(get(gcbo, ''string''));' footer ];
    cbroll    = [ header 'dattmp.transform(5) = str2num(get(gcbo, ''string''));' footer ];
    cbyaw     = [ header 'dattmp.transform(6) = str2num(get(gcbo, ''string''));' footer ];
    cbresizex = [ header 'dattmp.transform(7) = str2num(get(gcbo, ''string''));' footer ];
    cbresizey = [ header 'dattmp.transform(8) = str2num(get(gcbo, ''string''));' footer ];
    cbresizez = [ header 'dattmp.transform(9) = str2num(get(gcbo, ''string''));' footer ];
    cb_ok     = 'set(gcbo, ''userdata'', ''ok'')';
    cb_warp   = 'coregister(''warp'', gcbf);';
    cb_fid    = 'coregister(''fiducials'', gcbf);';
    
    opt = { 'unit', 'normalized', 'position' };
    h = uicontrol( opt{:}, [0    .15  1  .02], 'style', 'text', 'string', '');
    h = uicontrol( opt{:}, [0    .1  .2  .05], 'style', 'text', 'string', 'Move right {mm}');
    h = uicontrol( opt{:}, [0    .05 .2  .05], 'style', 'text', 'string', 'Move front {mm}' );
    h = uicontrol( opt{:}, [0    0   .2  .05], 'style', 'text', 'string', 'Move up {mm}');
    h = uicontrol( opt{:}, [0.2  .1  .1  .05], 'tag', 'right'  , 'callback', cbright  , 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.2  .05 .1  .05], 'tag', 'forward', 'callback', cbforward, 'style', 'edit', 'string', '' );
    h = uicontrol( opt{:}, [0.2  0   .1  .05], 'tag', 'up'     , 'callback', cbup     , 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.3  .1  .15 .05], 'style', 'text', 'string', 'Pitch (rad)');
    h = uicontrol( opt{:}, [0.3  .05 .15 .05], 'style', 'text', 'string', 'Roll (rad)' );
    h = uicontrol( opt{:}, [0.3  0   .15 .05], 'style', 'text', 'string', 'Yaw (rad)');
    h = uicontrol( opt{:}, [0.45 .1  .1  .05], 'tag', 'pitch', 'callback', cbpitch, 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.45 .05 .1  .05], 'tag', 'roll' , 'callback', cbroll , 'style', 'edit', 'string', '' );
    h = uicontrol( opt{:}, [0.45 0   .1  .05], 'tag', 'yaw'  , 'callback', cbyaw  , 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.55 .1  .15 .05], 'style', 'text', 'string', 'Resize {x}');
    h = uicontrol( opt{:}, [0.55 .05 .15 .05], 'style', 'text', 'string', 'Resize {y}' );
    h = uicontrol( opt{:}, [0.55 0   .15 .05], 'style', 'text', 'string', 'Resize {z}');
    h = uicontrol( opt{:}, [0.7  .1  .1  .05], 'tag', 'resizex', 'callback', cbresizex, 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.7  .05 .1  .05], 'tag', 'resizey', 'callback', cbresizey, 'style', 'edit', 'string', '' );
    h = uicontrol( opt{:}, [0.7  0   .1  .05], 'tag', 'resizez', 'callback', cbresizez, 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.8  .1  .2  .05], 'style', 'pushbutton', 'string', 'Align fiducials', 'callback', cb_fid);
    h = uicontrol( opt{:}, [0.8  .05 .2  .05], 'style', 'pushbutton', 'string', 'Warp montage', 'callback', cb_warp );
    h = uicontrol( opt{:}, [0.8  0   .1  .05], 'style', 'pushbutton', 'string', 'Cancel', 'callback', 'close(gcbf);' );
    h = uicontrol( opt{:}, [0.9  0   .1  .05], 'style', 'pushbutton', 'string', 'Ok', 'tag', 'ok', 'callback', cb_ok);
    
    % put labels next to electrodes
    % -----------------------------
    cb_label1 = [   'tmp = get(gcbf, ''userdata'');' ...
                    'if tmp.label1, set(gcbo, ''string'', ''Labels on'');' ...
                    'else           set(gcbo, ''string'', ''Labels off'');' ...
                    'end;' ...
                    'tmp.label1 = ~tmp.label1;' ...
                    'set(gcbf, ''userdata'', tmp);' ...
                    'clear tmp;' ...
                    'coregister(''redraw'', gcbf);' ];
    cb_label2 = [   'tmp = get(gcbf, ''userdata'');' ...
                    'if tmp.label2, set(gcbo, ''string'', ''Labels on'');' ...
                    'else           set(gcbo, ''string'', ''Labels off'');' ...
                    'end;' ...
                    'tmp.label2 = ~tmp.label2;' ...
                    'set(gcbf, ''userdata'', tmp);' ...
                    'clear tmp;' ...
                    'coregister(''redraw'', gcbf);' ];
    cb_mesh   = [   'tmp = get(gcbf, ''userdata'');' ...
                    'if tmp.meshon, set(gcbo, ''string'', ''Mesh on'');' ...
                    'else           set(gcbo, ''string'', ''Mesh off'');' ...
                    'end;' ...
                    'tmp.meshon = ~tmp.meshon;' ...
                    'set(gcbf, ''userdata'', tmp);' ...
                    'clear tmp;' ...
                    'coregister(''redraw'', gcbf);' ];
    cb_elecshow1 = [ 'tmp = get(gcbf, ''userdata'');' ...
                      'tmp.elecshow1 = pop_chansel( tmp.elec1.label, ''select'', tmp.elecshow1 );' ...
                      'if ~isempty(tmp.elecshow1), set(gcbf, ''userdata'', tmp);' ...
                      'coregister(''redraw'', gcbf); end; clear tmp;' ];
    cb_elecshow2 = [ 'tmp = get(gcbf, ''userdata'');' ...
                      'tmpstrs = { ''21 elec (10/20 system)'' ''86 elec (10/10 system)'' ''all elec (10/5 system)'' };' ...
                      'tmpres = inputgui( ''uilist'', {{ ''style'' ''text'' ''string'' ''show only'' } ' ...
                                        ' { ''style'' ''listbox'' ''string'' strvcat(tmpstrs) }}, ' ...
                                        ' ''geometry'', { 1 1 }, ''geomvert'', [1 3] );' ...
                      'if ~isempty(tmpres), tmp.elecshow2 = tmpstrs{tmpres{1}}; end;' ... 
                      'set(gcbf, ''userdata'', tmp);' ...
                      'clear tmp tmpres;' ...
                      'coregister(''redraw'', gcbf);' ];
    h = uicontrol( opt{:}, [0 0.75  .13 .05], 'style', 'pushbutton', 'string', 'Mesh off', 'callback', cb_mesh );

    % help message
    % ------------
    cb_helpme = [ 'warndlg2(strvcat( ''User channels (sometimes hidden in the 3-D mesh) are in green, reference channels in brown.'',' ...
            '''Press "Warp" to automatically warp user channels to corresponding reference channels.'',' ...
            '''Then, if desired, further edit the transformation manually using the coregister gui buttons.'',' ...
            ''' '',' ...
            '''To use locations of corresponding reference channels (and discard current locations),'',' ...
            '''select "Edit > Channel locations" in the EEGLAB mensu and press, "Look up loc." Select a '',' ...
            '''head model. Then re-open "Tools > Locate dipoles using DIPFIT2 > Head model and settings"'',' ...
            '''in the EEGLAB menu and select the "No coreg" option.'',' ];
    if ~isstruct(chanlocs2)
        if ~isempty(findstr(lower(chanlocs2), 'standard-10-5-cap385')) || ...
                ~isempty(findstr(lower(chanlocs2), 'standard_1005'))
            cb_helpme = [ cb_helpme '''Then re-open "Tools > Locate dipoles using DIPFIT2 > Head model and settings"'',' ...
                          '''in the EEGLAB menu and select the "No coreg" option.''), ''Warning'');' ];
        else
            cb_helpme = [ cb_helpme '''Then re-open the graphic interface function you were using.''), ''Warning'');' ];
        end
    end
    h = uicontrol( opt{:}, [0.87 0.95 .13 .05], 'style', 'pushbutton', 'string', 'Help me', 'callback',  cb_helpme);
    h = uicontrol( opt{:}, [0.87 0.90 .13 .05], 'style', 'pushbutton', 'string', 'Funct. help', 'callback', 'pophelp(''coregister'');' );

    % change colors
    % -------------
    hh = findobj('parent', gcf, 'style', 'text');
    set(hh, 'Backgroundcolor', GUIBACKCOLOR);
    set(hh, 'foregroundcolor', GUITEXTCOLOR);
    hh = findobj('parent', gcf, 'style', 'edit');
    set(hh, 'Backgroundcolor', GUIBACKCOLOR);
    set(hh, 'foregroundcolor', GUITEXTCOLOR);
    hh = findobj('parent', gcf, 'style', 'pushbutton');
    set(hh, 'Backgroundcolor', GUIBACKCOLOR);
    set(hh, 'foregroundcolor', GUITEXTCOLOR);

    h = uicontrol( opt{:}, [0  0.95 .13 .05], 'style', 'pushbutton', 'backgroundcolor', dat.color1, 'string', 'Labels on', 'callback', cb_label1 );
    h = uicontrol( opt{:}, [0  0.9  .13 .05], 'style', 'pushbutton', 'backgroundcolor', dat.color1, 'string', 'Electrodes', 'callback', cb_elecshow1 );

    h = uicontrol( opt{:}, [0 0.85  .13 .05], 'style', 'pushbutton', 'backgroundcolor', dat.color2, 'string', 'Labels on', 'callback', cb_label2 );
    h = uicontrol( opt{:}, [0 0.8   .13 .05], 'style', 'pushbutton', 'backgroundcolor', dat.color2, 'string', 'Electrodes', 'callback', cb_elecshow2 );
        
end

coregister('redraw', fid);
try icadefs; set(gcf, 'color', BACKCOLOR); catch, end

% wait until button press and return values
% -----------------------------------------
waitfor( findobj('parent', fid, 'tag', 'ok'), 'userdata');
try
    tmpobj = findobj(fid); % figure still exist ?
    if isempty(tmpobj), error(' '); end % After MATLAB 2014b
catch, transformmat = []; chanlocs1 = []; return; end
dat = get(fid, 'userdata');
transformmat = dat.transform;
chanlocs1    = dat.electransf;
title(g.title);
if ~strcmpi(g.manual, 'show')
    close(fid);
end

% plot electrodes
% ---------------
function plotelec(elec, elecshow, color, tag)
    
    X1 = elec.pnt(elecshow,1);
    Y1 = elec.pnt(elecshow,2);
    Z1 = elec.pnt(elecshow,3);    

    XL = xlim;
    YL = ylim;
    ZL = zlim;
    
    lim=max(1.05*max([X1;Y1;Z1])); %, max([XL YL ZL]));
    eps=lim/20;
    delete(findobj(gcf, 'tag', tag));
    
    % make bigger if fiducial
    % ------------------------
    fidlist = { 'nz' 'lpa' 'rpa' 'nazion' 'left' 'right' 'nasion' 'fidnz' 'fidt9' 'fidt10'};
    [~, fids ] = intersect_bc(lower(elec.label(elecshow)), fidlist);
    nonfids     = setdiff_bc(1:length(elec.label(elecshow)), fids);
    h1 = plot3(X1(nonfids),Y1(nonfids),Z1(nonfids), 'o', 'color', color); hold on;
    set(h1, 'tag', tag, 'marker', '.', 'markersize', 20);
    if ~isempty(fids)
        h2 = plot3(X1(fids),Y1(fids),Z1(fids), 'o', 'color', color*2/3); hold on;
        set(h2, 'tag', tag, 'marker', '.', 'markersize', 35); % make bigger if fiducial
    end

    % plot axis and labels
    %- -------------------
    if isempty(findobj(gcf, 'tag', 'axlabels'))
        plot3([0.08 0.12],[0 0],[0 0],'r','LineWidth',4) % nose
        plot3([0 lim],[0 0],[0 0],'b--', 'tag', 'axlabels')                 % axis
        plot3([0 0],[0 lim],[0 0],'g--')
        plot3([0 0],[0 0],[0 lim],'r--')
        plot3(0,0,0,'b+')
        text(lim+eps,0,0,'X','HorizontalAlignment','center',...
             'VerticalAlignment','middle','Color',[0 0 0],...
             'FontSize',10)
        text(0,lim+eps,0,'Y','HorizontalAlignment','center',...
             'VerticalAlignment','middle','Color',[0 0 0],...
             'FontSize',10)
        text(0,0,lim+eps,'Z','HorizontalAlignment','center',...
             'VerticalAlignment','middle','Color',[0 0 0],...
             'FontSize',10)
        box on
    end
    lim = abs(lim(1)); axis([-lim lim -lim lim -lim*0.5 lim]);
    axis equal;

% decode labels for electrode caps
% --------------------------------
% 86 channels
function indices = decodelabels( chanlocs, strchan );
    label1020 = { 'nz' 'lpa' 'rpa' 'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8',  'T7', 'C3', 'Cz', 'C4', 'T8',  'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'Oz', 'O2'}'; % 21 channels
    label1010 = { 'nz' 'lpa' 'rpa' 'Fp1', 'Fpz', 'Fp2', 'AF9', 'AF7', 'AF5', 'AF3', 'AF1', 'AFz', 'AF2', 'AF4', 'AF6', 'AF8', 'AF10', 'F9', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8', 'F10', 'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', 'T9', 'T7', 'C5', 'C3', 'C1', 'Cz', 'C2', ...
             'C4', 'C6', 'T8', 'T10', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10', 'P9', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO9', 'PO7', 'PO5', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO6', 'PO8', 'PO10', 'O1', 'Oz', 'O2', 'I1', 'Iz', 'I2'}'; ...
    if ~ischar(strchan), indices = strchan; return; end
    switch strchan
     case '21 elec (10/20 system)', indices = pop_chancoresp( struct('labels', chanlocs.label), struct('labels', label1020), 'gui', 'off');
     case '86 elec (10/10 system)', indices = pop_chancoresp( struct('labels', chanlocs.label), struct('labels', label1010), 'gui', 'off');
     case 'all elec (10/5 system)', indices = 1:length(chanlocs.label);
    otherwise, error('Unknown option');
  end
    
% plot electrode labels
% ---------------------
function plotlabels(elec, elecshow, color, tag);

    for i = 1:length(elecshow)
        coords = elec.pnt(elecshow(i),:);
        coords = coords*1.07;
        text(coords(1), coords(2), coords(3), elec.label{elecshow(i)},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',color, 'FontSize',10, 'tag', tag)
    end
    
% align fiducials
% ---------------
function [elec1, transf] = align_fiducials(elec1, elec2, fidnames1, fidnames2)

    % rename fiducials
    % ----------------
    ind1 = strmatch(fidnames1{1}, elec1.label, 'exact'); elec1.label{ind1} = fidnames2{1};
    ind2 = strmatch(fidnames1{2}, elec1.label, 'exact'); elec1.label{ind2} = fidnames2{2};
    ind3 = strmatch(fidnames1{3}, elec1.label, 'exact'); elec1.label{ind3} = fidnames2{3};
    cfg          = [];
    cfg.elec     = elec1;
    cfg.template = elec2;
    cfg.method   = 'realignfiducial'; 
    cfg.fiducial = fidnames2;
    elec3 = electroderealign(cfg);
    transf = homogenous2traditional(elec3.m);
    
    % test difference
    % ---------------
    diff1 = mean(mean(abs(elec3.m-traditionaldipfit(transf))));
    transf(6) = -transf(6);
    diff2 = mean(mean(abs(elec3.m-traditionaldipfit(transf))));
    if diff1 < diff2, transf(6) = -transf(6); end
    
    diff1 = mean(mean(abs(elec3.m-traditionaldipfit(transf))));
    transf(5) = -transf(5);
    diff2 = mean(mean(abs(elec3.m-traditionaldipfit(transf))));
    if diff1 < diff2, transf(5) = -transf(5); end
    
    diff1 = mean(mean(abs(elec3.m-traditionaldipfit(transf))));
    transf(4) = -transf(4);
    diff2 = mean(mean(abs(elec3.m-traditionaldipfit(transf))));
    if diff1 < diff2, transf(4) = -transf(4); end

    % rescale if necessary
    % --------------------
    coords1 = elec1.pnt([ind1 ind2 ind3],:); dist_coords1 = sqrt(sum(coords1.^2,2));
    ind1 = strmatch(fidnames2{1}, elec2.label, 'exact');
    ind2 = strmatch(fidnames2{2}, elec2.label, 'exact');
    ind3 = strmatch(fidnames2{3}, elec2.label, 'exact');
    coords2 = elec2.pnt([ind1 ind2 ind3],:); dist_coords2 = sqrt(sum(coords2.^2,2));
    ratio = mean(dist_coords2./dist_coords1);
    transf(7:9) = ratio;
    
    transfmat = traditionaldipfit(transf);
    elec1.pnt = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
    elec1.pnt = elec1.pnt(1:3,:)';
    
% warp channels
% -------------
function [elec1, transf] = warp_chans(elec1, elec2, chanlist, warpmethod)
    cfg          = [];
    cfg.elec     = elec1;
    cfg.template = elec2;
    cfg.method   = warpmethod;
    %cfg.feedback = 'yes';
    cfg.channel  = chanlist;
    if length(cfg.channel) > length(elec1.label)
        [ind] = intersect_bc(elec1.label, cfg.channel);
        elec1tmp.label   = elec1.label(ind);
        elec1tmp.elecpos = elec1.elecpos(ind,:);
        elec1tmp.pnt     = elec1.pnt(ind,:);
        cfg.elec         = elec1tmp;
    end
    try
        elec3 = electroderealign(cfg);
    catch
        error( [ 'You need to select more pairs or the correspondance you selectd' 10 ...
                 'leads to a failure in initial objective function evaluation' 10 ...
                 'which means that the correspondance is wrong.' ]);
    end
    [tmp ind1 ] = intersect_bc( lower(elec1.label), lower(chanlist) );
    [tmp ind2 ] = intersect_bc( lower(elec2.label), lower(chanlist) );
    
    transf = elec3.m;
    transf(4:6) = transf(4:6)/180*pi;
    if length(transf) == 6, transf(7:9) = 1; end
    transf = checktransf(transf, elec1, elec2);
    
    dpre = mean(sqrt(sum((elec1.pnt(ind1,:) - elec2.pnt(ind2,:)).^2, 2)));
    transfmat = traditionaldipfit(transf);
    elec1.pnt = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
    elec1.pnt = elec1.pnt(1:3,:)';
    dpost = mean(sqrt(sum((elec1.pnt(ind1,:) - elec2.pnt(ind2,:)).^2, 2)));
    
    fprintf('mean distance prior to warping %f, after warping %f\n', dpre, dpost);
    
% test difference and invert axis if necessary
% --------------------------------------------
function transf = checktransf(transf, elec1, elec2)
    
    [tmp ind1 ind2] = intersect_bc( elec1.label, elec2.label );
    
    transfmat = traditionaldipfit(transf);
    tmppnt = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
    tmppnt = tmppnt(1:3,:)';
    diff1  = tmppnt(ind1,:) - elec2.pnt(ind2,:);
    diff1  = mean(sum(diff1.^2,2));
    
    transf(6) = -transf(6); % yaw angle is sometimes inverted
    transfmat = traditionaldipfit(transf);
    tmppnt = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
    tmppnt = tmppnt(1:3,:)';
    diff2  = tmppnt(ind1,:) - elec2.pnt(ind2,:);
    diff2  = mean(sum(diff2.^2,2));
    
    if diff1 < diff2, transf(6) = -transf(6); else diff1 = diff2; end

    transf(4) = -transf(4); % yaw angle is sometimes inverted
    transfmat = traditionaldipfit(transf);
    tmppnt = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
    tmppnt = tmppnt(1:3,:)';
    diff2  = tmppnt(ind1,:) - elec2.pnt(ind2,:);
    diff2  = mean(sum(diff2.^2,2));
    
    if diff1 < diff2, transf(4) = -transf(4); end
    
% redraw GUI
% ----------
function redrawgui(fid)
    dat = get(fid, 'userdata');
    tmpobj = findobj(fid, 'tag', 'pitch'); set(tmpobj, 'string', num2str(dat.transform(4),4));
    tmpobj = findobj(fid, 'tag', 'roll' ); set(tmpobj, 'string', num2str(dat.transform(5),4));
    tmpobj = findobj(fid, 'tag', 'yaw'  ); set(tmpobj, 'string', num2str(dat.transform(6),4));
    tmpobj = findobj(fid, 'tag', 'right'  ); set(tmpobj, 'string', num2str(dat.transform(1),4));
    tmpobj = findobj(fid, 'tag', 'forward'); set(tmpobj, 'string', num2str(dat.transform(2),4));
    tmpobj = findobj(fid, 'tag', 'up'     ); set(tmpobj, 'string', num2str(dat.transform(3),4));
    tmpobj = findobj(fid, 'tag', 'resizex'); set(tmpobj, 'string', num2str(dat.transform(7),4));
    tmpobj = findobj(fid, 'tag', 'resizey'); set(tmpobj, 'string', num2str(dat.transform(8),4));
    tmpobj = findobj(fid, 'tag', 'resizez'); set(tmpobj, 'string', num2str(dat.transform(9),4));
    tmpview = view;
    
    if size(dat.transform,1) > 1
        dat.electransf.pnt = dat.transform*[ dat.elec1.pnt ones(size(dat.elec1.pnt,1),1) ]';
    else
        dat.electransf.pnt = traditionaldipfit(dat.transform)*[ dat.elec1.pnt ones(size(dat.elec1.pnt,1),1) ]';
    end
    dat.electransf.pnt   = dat.electransf.pnt(1:3,:)';
    dat.electransf.label = dat.elec1.label;
    set(fid, 'userdata', dat);
    
    h = findobj(fid, 'tag', 'plot3d');
    if isempty(h)
        axis off;
        h = axes('unit', 'normalized', 'position', [0 0.2 1 0.75]);
        set(h, 'tag', 'plot3d');
        axis off;
    else 
        axes(h);
        %axis off;
    end
    plotelec(dat.electransf, dat.elecshow1, dat.color1, 'elec1');
    if ~isempty(dat.elec2)
        dat.elecshow2 = decodelabels( dat.elec2, dat.elecshow2 );
        plotelec(dat.elec2, dat.elecshow2, dat.color2, 'elec2');
    end
    set(h, 'tag', 'plot3d');
    
    % plot mesh
    % ---------
    if ~isempty(dat.meshpnt) && isempty(findobj(gcf, 'tag', 'mesh'))
        if ~isempty(dat.meshtri)
            p1 = plotmesh(dat.meshtri, dat.meshpnt, [], 1);
            set(p1, 'tag', 'mesh');
        else
            facecolor(1,1,1) = 1; facecolor(1,1,2) = .75; facecolor(1,1,3) = .65;
            cdat = repmat( facecolor, [ size(dat.meshpnt{1}) 1]);
            h = mesh(dat.meshpnt{1}, dat.meshpnt{2}, dat.meshpnt{3}, ...
                'cdata', cdat, 'tag', 'mesh', 'facecolor', squeeze(facecolor), 'edgecolor', 'none');
            hidden off;
            lightangle(45,30);
            lightangle(45+180,30);
            lighting phong
            s = plotnose([85 0 -75 0 0 pi/2 10 10 40]);
            set(s, 'tag', 'mesh');
        end
    end
    meshobj = findobj(gcf, 'tag', 'mesh');
    if dat.meshon
         set( meshobj, 'visible', 'on');
    else set( meshobj, 'visible', 'off');
    end
    
    % plot electrodes
    % ---------------
    delete(findobj(gcf, 'tag', 'elec1labels'));        
    delete(findobj(gcf, 'tag', 'elec2labels'));        
    if dat.label1
        plotlabels(dat.electransf, dat.elecshow1, dat.color1, 'elec1labels');
    end
    if dat.label2
        plotlabels(dat.elec2, dat.elecshow2, dat.color2*0.5, 'elec2labels');
    end
    
    %view(tmpview);
    rotate3d on    
    
  
% function to plot the nose
% -------------------------
function s = plotnose(transf, col)

    if nargin < 1
        transf = [0 0 0 0 0 0 1 1 1];
    end
    if nargin < 2
        col = [1 0.75 0.65 ];
    end
        
    x=[ % cube
     NaN -1 1 NaN
      -1 -1 1 1
      -1 -1 1 1
     NaN -1 1 NaN
     NaN -1 1 NaN
     NaN NaN NaN NaN
     ];
    
    y=[ % cube
     NaN -1 -1   NaN
      -1 -1 -1   -1
       1  1  1    1
     NaN  1  1   NaN
     NaN -1 -1   NaN
     NaN NaN NaN NaN
     ];
    
    z=[ % cube
     NaN 0 0 NaN
       0 1 1 0
       0 1 1 0
     NaN 0 0 NaN
     NaN 0 0 NaN
     NaN NaN NaN NaN
     ];

    x=[ % noze
     NaN -1  1 NaN
      -1  0  0 1
     -.3  0  0 .3
     NaN -.3 .3 NaN
     NaN -1  1 NaN
     NaN NaN NaN NaN
     ];
    
    y=[ % noze
     NaN -1 -1   NaN
      -1 -1 -1   -1
       1 -1 -1    1
     NaN  1  1   NaN
     NaN -1 -1   NaN
     NaN NaN NaN NaN
     ];
    
    z=[ % noze
     NaN 0 0 NaN
       0 1 1 0
       0 1 1 0
     NaN 0 0 NaN
     NaN 0 0 NaN
     NaN NaN NaN NaN
     ];

    % apply homogenous transformation
    % -------------------------------
    transfhom = traditionaldipfit( transf );
    xyz       = [ x(:) y(:) z(:) ones(length(x(:)),1) ];
    xyz2      = transfhom * xyz';
    x(:)      = xyz2(1,:)';
    y(:)      = xyz2(2,:)';
    z(:)      = xyz2(3,:)';

    % dealing with colors
    % -------------------
    cc=zeros(8,3);
    cc(1,:) = col;
    cc(2,:) = col;
    cc(3,:) = col;
    cc(4,:) = col;
    cc(5,:) = col;
    cc(6,:) = col;
    cc(7,:) = col;
    cc(8,:) = col;
    cc(1,:)=[0 0 0]; % black
    cc(2,:)=[1 0 0]; % red
    cc(3,:)=[0 1 0]; % green
    cc(4,:)=[0 0 1]; % blue
    cc(5,:)=[1 0 1]; % magenta
    cc(6,:)=[0 1 1]; % cyan
    cc(7,:)=[1 1 0]; % yellow
    cc(8,:)=[1 1 1]; % white
    cs=size(x);
    c=repmat(zeros(cs),[1 1 3]);
    for i=1:size(cc,1)
        ix=find(x==cc(i,1) &...
                y==cc(i,2) &...
                z==cc(i,3));
        [ir,ic]=ind2sub(cs,ix);
        for k=1:3
            for m=1:length(ir)
                c(ir(m),ic(m),k)=cc(i,k);
            end
        end
    end
    
    % plotting surface
    % ----------------
    facecolor = zeros(size(x,1), size(z,2), 3);
    facecolor(:,:,1) = 1; facecolor(:,:,2) = .75; facecolor(:,:,3) = .65;
   
    s=surf(x,y,z,facecolor);
    set(s, 'edgecolor', [0.5 0.5 0.5]);

