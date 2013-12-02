% pop_dipfit_settings() - select global settings for dipole fitting through a pop up window
%
% Usage:
%   >> OUTEEG = pop_dipfit_settings ( INEEG ); % pop up window
%   >> OUTEEG = pop_dipfit_settings ( INEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
% Inputs:
%   INEEG	input dataset
%
% Optional inputs:
%   'hdmfile'  - [string] file containing a head model compatible with
%                the Fieldtrip dipolefitting() function ("vol" entry)
%   'mrifile'  - [string] file containing an anatomical MR head image. 
%                The MRI must be normalized to the MNI brain. See the .mat 
%                files used by the sphere and boundary element models
%                (For instance, select the sphere model and study 'EEG.dipfit'). 
%                If SPM2 software is installed, dipfit will be able to read 
%                most MRI file formats for plotting purposes (.mnc files, etc...). 
%                To plot dipoles in a subject MRI, first normalize the MRI 
%                to the MNI brain using SPM2.
%   'coordformat' - ['MNI'|'Spherical'] Coordinates returned by the selected
%                head model. May be MNI coordinates or spherical coordinates
%                (For spherical coordinates, the head radius is assumed to be 85 mm.
%   'chanfile' - [string] template channel locations file. (This function will
%                check whether your channel locations file is compatible with 
%                your selected head model).
%   'chansel'  - [integer vector] indices of channels to use for dipole fitting. 
%                {default: all}
%   'coord_transform' - [float array] Talairach transformation matrix for
%                       aligning the dataset channel locations to the selected 
%                       head model.
%   'electrodes'      - [integer array] indices of channels to include
%                       in the dipole model. {default: all}
% Outputs:
%   OUTEEG	output dataset
%
% Author: Arnaud Delorme, SCCN, La Jolla 2003-
%         Robert Oostenveld, SMI/FCDC, Nijmegen 2003

% MEG flag:
%   'gradfile' - [string] file containing gradiometer locations
%                ("gradfile" parameter in Fieldtrip dipolefitting() function)

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl

% Copyright (C) 2003 arno@salk.edu, Arnaud Delorme, SCCN, La Jolla 2003-2005
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUTEEG, com] = pop_dipfit_settings ( EEG, varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help pop_dipfit_settings;
   return;
end;

if ~plugin_askinstall('Fieldtrip-lite', 'ft_sourceanalysis'), return; end;

OUTEEG = EEG;
com = '';

% get the default values and filenames
dipfitdefs;

if nargin < 2
    
    if isstr(EEG) % setmodel
        tmpdat = get(gcf, 'userdata');
        chanfile = tmpdat.chanfile;
        tmpdat   = tmpdat.template_models;
        tmpval = get(findobj(gcf, 'tag', 'listmodels'), 'value');
        set(findobj(gcf, 'tag', 'model'), 'string', char(tmpdat(tmpval).hdmfile));
        set(findobj(gcf, 'tag', 'coord'), 'value' , fastif(strcmpi(tmpdat(tmpval).coordformat,'MNI'),2, ...
                                                     fastif(strcmpi(tmpdat(tmpval).coordformat,'CTF'),3,1)));
        set(findobj(gcf, 'tag', 'mri'  ), 'string', char(tmpdat(tmpval).mrifile));
        set(findobj(gcf, 'tag', 'meg'), 'string', char(tmpdat(tmpval).chanfile));
        set(findobj(gcf, 'tag', 'coregcheckbox'), 'value', 0);
        if tmpval < 3,
            set(findobj(gcf, 'userdata', 'editable'), 'enable', 'off');
        else,
            set(findobj(gcf, 'userdata', 'editable'), 'enable', 'on');
        end;
        if tmpval == 3,
            set(findobj(gcf, 'tag', 'headstr'), 'string', 'Subject CTF head model file (default.htm)');
            set(findobj(gcf, 'tag', 'mristr'),  'string', 'Subject MRI (coregistered with CTF head)');
            set(findobj(gcf, 'tag', 'chanstr'), 'string', 'CTF Res4 file');
            set(findobj(gcf, 'tag', 'manualcoreg'), 'enable', 'off');
            set(findobj(gcf, 'userdata', 'coreg'), 'enable', 'off');
        else,
            set(findobj(gcf, 'tag', 'headstr'), 'string', 'Head model file');
            set(findobj(gcf, 'tag', 'mristr'),  'string', 'MRI file');
            set(findobj(gcf, 'tag', 'chanstr'), 'string', 'Model template channel locations file');
            set(findobj(gcf, 'tag', 'manualcoreg'), 'enable', 'on');
            set(findobj(gcf, 'userdata', 'coreg'), 'enable', 'on');
        end;
        tmpl = tmpdat(tmpval).coord_transform;
        set(findobj(gcf, 'tag', 'coregtext'), 'string', '');
        set(findobj(gcf, 'tag', 'coregcheckbox'), 'value', 0);
        [allkeywordstrue transform] = lookupchantemplate(chanfile, tmpl);
        if allkeywordstrue,
            set(findobj(gcf, 'tag', 'coregtext'), 'string', char(vararg2str({ transform })));
            if isempty(transform)
                 set(findobj(gcf, 'tag', 'coregcheckbox'), 'value', 1);
            else set(findobj(gcf, 'tag', 'coregcheckbox'), 'value', 0);
            end;
        end;
        return;
    end;
    
    % detect DIPFIT1.0x structure
    % ---------------------------
    if isfield(EEG.dipfit, 'vol')
        str = [ 'Dipole information structure from DIPFIT v1.02 detected.' ...
                'Keep or erase the old dipole information including dipole locations? ' ...
                'In either case, a new dipole model can be constructed.' ];
        
        tmpButtonName=questdlg2( strmultiline(str, 60), 'Old DIPFIT structure', 'Keep', 'Erase', 'Keep');
        if strcmpi(tmpButtonName, 'Keep'), return; end;       

    elseif isfield(EEG.dipfit, 'hdmfile')
        % detect previous DIPFIT structure
        % --------------------------------
        str = [ 'Dipole information and settings are present in the dataset. ' ...
                'Keep or erase this information?' ];
        tmpButtonName=questdlg2( strmultiline(str, 60), 'Old DIPFIT structure', 'Keep', 'Erase', 'Keep');
        if strcmpi(tmpButtonName, 'Keep'), return; end;       
    end;    
    
    % define the callbacks for the buttons
    % -------------------------------------
    cb_selectelectrodes = [ 'tmplocs = EEG.chanlocs; tmp = select_channel_list({tmplocs.label}, ' ...
                            'eval(get(findobj(gcbf, ''tag'', ''elec''), ''string'')));' ...
                            'set(findobj(gcbf, ''tag'', ''elec''), ''string'',[''[''  num2str(tmp) '']'']); clear tmplocs;' ]; % did not work
    cb_selectelectrodes = 'tmplocs = EEG.chanlocs; set(findobj(gcbf, ''tag'', ''elec''), ''string'', int2str(pop_chansel({tmplocs.labels}))); clear tmplocs;';
    cb_volmodel = [ 'tmpdat = get(gcbf, ''userdata'');' ... 
                    'tmpind = get(gcbo, ''value'');' ... 
                    'set(findobj(gcbf, ''tag'', ''radii''),   ''string'', num2str(tmpdat{tmpind}.r,3));' ...
                    'set(findobj(gcbf, ''tag'', ''conduct''), ''string'', num2str(tmpdat{tmpind}.c,3));' ...
                    'clear tmpdat tmpind;' ];
    cb_changeradii   = [  'tmpdat = get(gcbf, ''userdata'');' ...
                          'tmpdat.vol.r = str2num(get(gcbo, ''string''));' ...
                          'set(gcf, ''userdata'', tmpdat)' ];
    cb_changeconduct = [  'tmpdat = get(gcbf, ''userdata'');' ...
                          'tmpdat.vol.c = str2num(get(gcbo, ''string''));' ...
                          'set(gcf, ''userdata'', tmpdat)' ];
    cb_changeorigin  = [  'tmpdat = get(gcbf, ''userdata'');' ...
                          'tmpdat.vol.o = str2num(get(gcbo, ''string''));' ...
                          'set(gcf, ''userdata'', tmpdat)' ];
    % cb_fitelec = [ 'if get(gcbo, ''value''),' ...
    %                '  set(findobj(gcbf, ''tag'', ''origin''), ''enable'', ''off'');' ...
    %                'else' ...
    %                '  set(findobj(gcbf, ''tag'', ''origin''), ''enable'', ''on'');' ...
    %                'end;' ];
    valmodel     = 1;
    userdata    = [];
    if isfield(EEG.chaninfo, 'filename')
        if ~isempty(findstr(lower(EEG.chaninfo.filename), 'standard-10-5-cap385')), valmodel = 1; end;
        if ~isempty(findstr(lower(EEG.chaninfo.filename), 'standard_1005')),        valmodel = 2; end;
    end;
   
    geomvert = [3 1 1 1 1 1 1 1 1 1 1];
    
    geomhorz = {
        [1 2] 
        [1]
        [1 1.3 0.5 0.5 ]
        [1 1.3 0.9 0.1 ]
        [1 1.3 0.5 0.5 ]
        [1 1.3 0.5 0.5 ]
        [1 1.3 0.5 0.5 ]
        [1 1.3 0.5 0.5 ]
        [1]
        [1]
        [1] };
    
    % define each individual graphical user element
    comhelp1 = [ 'warndlg2(strvcat(''The two default head models are in ''standard_BEM'' and ''standard_BESA'''',' ...
                 ''' sub-folders in the DIPFIT2 plugin folder, and may be modified there.''), ''Model type'');' ];
    comhelp3 = [ 'warndlg2(strvcat(''Any MR image normalized to the MNI brain model may be used for plotting'',' ...
                 '''(see the DIPFIT 2.0 tutorial for more information)''), ''Model type'');' ];
    comhelp2 = [ 'warndlg2(strvcat(''The template location file associated with the head model'',' ...
                 '''you are using must be entered (see tutorial).''), ''Template location file'');' ];
    commandload1 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''model''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    commandload2 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''meg''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    commandload3 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''mri''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    cb_selectcoreg = [ 'tmpmodel = get( findobj(gcbf, ''tag'', ''model''), ''string'');' ...
                       'tmploc2  = get( findobj(gcbf, ''tag'', ''meg'')  , ''string'');' ...
                       'tmploc1  = get( gcbo, ''userdata'');' ...
                       'tmptransf = get( findobj(gcbf, ''tag'', ''coregtext''), ''string'');' ...
                       '[tmp tmptransf] = coregister(tmploc1{1}, tmploc2, ''mesh'', tmpmodel,' ...
                       '                       ''transform'', str2num(tmptransf), ''chaninfo1'', tmploc1{2}, ''helpmsg'', ''on'');' ...
                       'if ~isempty(tmptransf), set( findobj(gcbf, ''tag'', ''coregtext''), ''string'', num2str(tmptransf)); end;' ...
                       'clear tmpmodel tmploc2 tmploc1 tmp tmptransf;' ];
    setmodel = [ 'pop_dipfit_settings(''setmodel'');' ];
    
    dipfitdefs; % contains template_model
        
    templatenames = { template_models.name };
    elements  = { ...
        { 'style' 'text'        'string'  [ 'Head model (click to select)' 10 '' ] } ...
        { 'style' 'listbox'     'string'  strvcat(templatenames{:}) ... 
                                'callback' setmodel 'value' valmodel 'tag' 'listmodels' } { } ...
        { 'style' 'text'        'string' '________' 'tag' 'headstr' } ...
        { 'style' 'edit'        'string' '' 'tag'      'model' 'userdata' 'editable' 'enable' 'off'} ...
        { 'style' 'pushbutton'  'string' 'Browse'    'callback' commandload1       'userdata' 'editable' 'enable' 'off' } ...
        { 'style' 'pushbutton'  'string' 'Help'      'callback' comhelp1 } ...
        { 'style' 'text'        'string' 'Output coordinates' } ...
        { 'style' 'popupmenu'   'string' 'spherical (head radius 85 mm)|MNI|CTF' 'tag' 'coord' ...
          'value' 1  'userdata' 'editable' 'enable' 'off'} ...
        { 'style' 'text'        'string' 'Click to select' } { } ...
        { 'style' 'text'        'string' '________' 'tag' 'mristr' } ...
        { 'style' 'edit'        'string' '' 'tag' 'mri' } ...
        { 'style' 'pushbutton'  'string' 'Browse'       'callback' commandload3 } ...
        { 'style' 'pushbutton'  'string' 'Help'         'callback' comhelp3 } ...
        { 'style' 'text'        'string' '________', 'tag', 'chanstr' } ...
        { 'style' 'edit'        'string' '' 'tag' 'meg'  'userdata' 'editable' 'enable' 'off'} ...
        { 'style' 'pushbutton'  'string' 'Browse'       'callback' commandload2  'userdata' 'editable' 'enable' 'off'} ...
        { 'style' 'pushbutton'  'string' 'Help'         'callback' comhelp2 } ...
        { 'style' 'text'        'string' 'Co-register chan. locs. with head model' 'userdata' 'coreg' } ...
        { 'style' 'edit'        'string' '' 'tag' 'coregtext' 'userdata' 'coreg' } ...
        { 'style' 'pushbutton'  'string' 'Manual Co-Reg.' 'tag' 'manualcoreg' 'callback' cb_selectcoreg 'userdata' { EEG.chanlocs,EEG.chaninfo } } ... 
        { 'style' 'checkbox'    'string' 'No Co-Reg.'    'tag' 'coregcheckbox' 'value' 0  'userdata' 'coreg' } ... 
        { 'style' 'text'        'string' 'Channels to omit from dipole fitting' } ...
        { 'style' 'edit'        'string' ''             'tag' 'elec' } ...
        { 'style' 'pushbutton'  'string' 'List' 'callback' cb_selectelectrodes } { } ...
        { } ...
        { 'style' 'text'        'string' 'Note: For EEG, check that the channel locations are on the surface of the head model' } ...
        { 'style' 'text'        'string' '(To do this: ''Set head radius'' to about 85 in the channel editor).' } ...
                };
    
    % plot GUI and protect parameters
    % -------------------------------
    userdata.template_models  = template_models;
    if isfield(EEG.chaninfo, 'filename')
         userdata.chanfile         = lower(EEG.chaninfo.filename);
    else userdata.chanfile         = '';
    end;
    optiongui = { 'geometry', geomhorz, 'uilist', elements, 'helpcom', 'pophelp(''pop_dipfit_settings'')', ...
                  'title', 'Dipole fit settings - pop_dipfit_settings()', ...
                  'userdata', userdata, 'geomvert', geomvert 'eval' 'pop_dipfit_settings(''setmodel'');' };
	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', 'noclose', optiongui{:});
    if isempty(result), return; end;
    if ~isempty(get(0, 'currentfigure')) currentfig = gcf; else return; end;
    
    while test_wrong_parameters(currentfig)
    	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', currentfig, optiongui{:});
        if isempty(result), return; end;
    end;
    close(currentfig);

    % decode GUI inputs
    % -----------------
    options = {};
    options = { options{:} 'hdmfile'      result{2} };
    options = { options{:} 'coordformat'  fastif(result{3} == 2, 'MNI', fastif(result{3} == 1, 'Spherical', 'CTF')) };
    options = { options{:} 'mrifile'      result{4} };
    options = { options{:} 'chanfile'     result{5} };
    if ~result{7}, options = { options{:} 'coord_transform' str2num(result{6}) }; end;
    options = { options{:} 'chansel'      setdiff(1:EEG.nbchan, str2num(result{8})) };

else
    options = varargin;
end

g = finputcheck(options, { 'hdmfile'  'string'    []         '';
                                 'mrifile'  'string'    []         '';
                                 'chanfile' 'string'    []         '';
                                 'chansel'  'integer'   []         [1:EEG.nbchan];
                                 'electrodes' 'integer'   []         [];
                                 'coord_transform' 'real' []         [];
                                 'coordformat' 'string'    { 'MNI','spherical','CTF' } 'MNI' });
if isstr(g), error(g); end;

OUTEEG = rmfield(OUTEEG, 'dipfit');
OUTEEG.dipfit.hdmfile     = g.hdmfile;
OUTEEG.dipfit.mrifile     = g.mrifile;
OUTEEG.dipfit.chanfile    = g.chanfile;
OUTEEG.dipfit.chansel     = g.chansel;
OUTEEG.dipfit.coordformat     = g.coordformat;
OUTEEG.dipfit.coord_transform = g.coord_transform;
if ~isempty(g.electrodes), OUTEEG.dipfit.chansel = g.electrodes; end;

% removing channels with no coordinates
% -------------------------------------
[tmpeloc labels Th Rd indices] = readlocs(EEG.chanlocs);
if length(indices) < length(EEG.chanlocs)
    disp('Warning: Channels removed from dipole fitting no longer have location coordinates!');
    OUTEEG.dipfit.chansel = intersect( OUTEEG.dipfit.chansel, indices);
end;

% checking electrode configuration
% --------------------------------
if 0
    disp('Checking the electrode configuration');
    tmpchan          = readlocs(OUTEEG.dipfit.chanfile);
    [tmp1 ind1 ind2] = intersect( lower({ tmpchan.labels }), lower({ OUTEEG.chanlocs.labels }));
    if isempty(tmp1)
        disp('No channel labels in common found between template and dataset channels');
        if ~isempty(findstr(OUTEEG.dipfit.hdmfile, 'BESA'))
            disp('Use the channel editor to fit a head sphere to your channel locations.');
            disp('Check for inconsistency in dipole info.');
        else
            disp('Results using standard BEM model are INACCURATE when the chan locations are not on the head surface!');
        end;
    else % common channels: performing best transformation
        TMP = OUTEEG;
        elec1 = eeglab2fieldtrip(TMP, 'elec');
        elec1 = elec1.elec;
        TMP.chanlocs = tmpchan;
        elec2 = eeglab2fieldtrip(TMP, 'elec');
        elec2 = elec2.elec;
        cfg.elec     = elec1;
        cfg.template = elec2;
        cfg.method   = 'warp';
        elec3 = electrodenormalize(cfg);
    
        % convert back to EEGLAB format
        OUTEEG.chanlocs = struct( 'labels', elec3.label, ...
                              'X'     , mat2cell(elec3.pnt(:,1)'), ...
                              'Y'     , mat2cell(elec3.pnt(:,2)'), ...
                              'Z'     , mat2cell(elec3.pnt(:,3)') );
        OUTEEG.chanlocs = convertlocs(OUTEEG.chanlocs, 'cart2all');
    end;
    
end;

com = sprintf('%s = pop_dipfit_settings( %s, %s);', inputname(1), inputname(1), vararg2str(options));

% test for wrong parameters
% -------------------------
function bool = test_wrong_parameters(hdl)

    coreg1 = get( findobj( hdl, 'tag', 'coregtext')    , 'string' );
    coreg2 = get( findobj( hdl, 'tag', 'coregcheckbox'), 'value' );
    meg    = get( findobj( hdl, 'tag', 'coord'), 'value' );
    
    bool = 0;
    if meg == 3, return; end;
    if coreg2 == 0 & isempty(coreg1)
         bool = 1; warndlg2(strvcat('You must co-register your channel locations', ...
                                    'with the head model (Press buttun, "Manual Co-Reg".', ...
                                    'and follow instructions); To bypass co-registration,', ...
                                    'check the checkbox " No Co-Reg".'), 'Error');
    end;
